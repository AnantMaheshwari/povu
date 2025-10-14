#include "povu/genomics/rov.hpp"
#include <algorithm> // for min, max
#include <atomic>    // for atomic, memory_order
#include <cmath>     // for ceil
#include <cstddef>   // for size_t
#include <cstdlib>
#include <iterator> // for back_insert_iterator, bac...
#include <map>	    // for map
#include <optional> // for optional, operator==
#include <string>   // for basic_string, string
#include <utility>  // for move
#include <vector>   // for vector

#include "fmt/core.h"		  // for format_to
#include "indicators/setting.hpp" // for PostfixText
#include "povu/common/app.hpp"	  // for config
#include "povu/common/constants.hpp"
#include "povu/common/core.hpp" // for pt
#include "povu/common/log.hpp"
#include "povu/common/progress.hpp" // for set_progress_bar_com...
#include "povu/common/utils.hpp"
#include "povu/genomics/graph.hpp" // for RoV, find_walks, pgt
// #include "povu/genomics/vcf.hpp"
#include "povu/graph/pvst.hpp" // for Tree, VertexBase
#include "povu/graph/types.hpp"

#include "povu/common/utils.hpp" // for print_with_comma

namespace povu::genomics::rov
{
namespace pgg = povu::genomics::graph;
namespace pvst = povu::pvst;
using namespace povu::progress;

constexpr var_type_e ins = var_type_e::ins;
constexpr var_type_e del = var_type_e::del;
constexpr var_type_e sub = var_type_e::sub;

// ------------
bool is_valid(pt::u32 i)
{
	return i != pc::INVALID_IDX;
}

/**
 * map each step to its position in the walk
 */
std::map<ptg::step_t, pt::u32> positions(const ptg::walk_t &w)
{
	std::map<ptg::step_t, pt::u32> m;
	for (pt::u32 i{}; i < w.size(); i++)
		m[w[i]] = i;

	return m;
}

/**
 * look up each step in the walk in the map
 */
std::vector<pt::u32> comp_lookup(const ptg::walk_t &w,
				 std::map<ptg::step_t, pt::u32> &m)
{
	std::vector<pt::u32> lu(w.size(), pc::INVALID_IDX);
	for (pt::u32 i{}; i < w.size(); i++)
		if (pv_cmp::contains(m, w[i]))
			lu[i] = m[w[i]];

	return lu;
}

pt::slice_t find_context(std::vector<pt::u32> &w, pt::u32 i)
{
	pt::u32 left = i;
	pt::u32 right = i + 1;

	while (right < w.size() && !is_valid(w[right]))
		right++;

	return {left, right - left};
}

std::pair<std::vector<pt::slice_t>, std::vector<var_type_e>>
find_rovs(std::vector<pt::u32> lu)
{
	std::vector<pt::slice_t> slices;
	std::vector<var_type_e> var_types;

	// because we move left to right we will always
	// start with the leftmost side of the slice
	for (pt::u32 i{}; i < lu.size();) {
		if (!is_valid(lu[i])) {
			pt::slice_t sl = find_context(lu, i);
			auto t = (i > 0) && (lu[i + sl.len] - lu[i - 1] == 1)
					 ? ins
					 : sub;

			var_types.push_back(t);
			slices.push_back(sl);
			i += sl.len;
			continue;
		}

		// del
		if (i > 0 && is_valid(lu[i - 1]) && is_valid(lu[i]) &&
		    lu[i] - lu[i - 1] != 1) {
			slices.emplace_back(i - 1, lu[i] - lu[i - 1] - 1);
			var_types.emplace_back(var_type_e::del);
		}

		i++;
	}

	return {std::move(slices), std::move(var_types)};
}

pt::op_t<std::pair<std::vector<pt::slice_t>, std::vector<var_type_e>>>
lineup_pairs(const ptg::walk_t &w1, const ptg::walk_t &w2)
{
	std::map<ptg::step_t, pt::u32> pos_map1 = positions(w1);
	std::map<ptg::step_t, pt::u32> pos_map2 = positions(w2);

	std::vector<pt::u32> lu1 = comp_lookup(w1, pos_map2);
	std::vector<pt::u32> lu2 = comp_lookup(w2, pos_map1);

	auto rovs1 = find_rovs(lu1);
	auto rovs2 = find_rovs(lu2);

	return {rovs1, rovs2};
}

void find_hidden(RoV &r)
{
	auto compare_pair = [&](pt::u32 i, pt::u32 j)
	{
		auto [a, b] = lineup_pairs(r.get_walk(i), r.get_walk(j));
		if (!a.first.empty())
			r.add_extra({{i, j}, a.first, a.second});
		if (!b.first.empty())
			r.add_extra({{j, i}, b.first, b.second});
	};

	std::set<pt::up_t<pt::u32>> seen;
	for (pt::u32 i{}; i < r.walk_count(); i++) {
		for (pt::u32 j{}; j < r.walk_count(); j++) {
			auto p = pt::up_t<pt::u32>{i, j};

			if (i == j || pv_cmp::contains(seen, p))
				continue; // skip self or already seen

			compare_pair(i, j);
			seen.insert(p);
		}
	}
}

// ------------

/**
 * Check if a vertex in the pvst is a flubble leaf
 * A flubble leaf is a vertex that has no children that are also
 * flubbles
 */
bool is_fl_leaf(const pvst::Tree &pvst, pt::idx_t pvst_v_idx) noexcept
{
	const pvst::VertexBase *pvst_v_ptr =
		pvst.get_vertex_const_ptr(pvst_v_idx);

	// we assume that the vertex has a clan
	pvst::vf_e prt_fam = pvst_v_ptr->get_fam();
	if (pvst::to_clan(prt_fam).value() != pvst::vc_e::fl_like) {
		return false; // not a flubble
	}

	for (pt::idx_t v_idx : pvst.get_children(pvst_v_idx)) {
		pvst::vf_e c_fam = pvst.get_vertex_const_ptr(v_idx)->get_fam();
		if (pvst::to_clan(c_fam) == pvst::vc_e::fl_like) {
			return false;
		}
	}

	return true;
}

/**
 * true when the vertex is a flubble leaf or a leaf in the pvst
 */
bool should_call(const pvst::Tree &pvst, const pvst::VertexBase *pvst_v_ptr,
		 pt::idx_t pvst_v_idx)
{
	if (auto opt_route_params = pvst_v_ptr->get_route_params())
		return is_fl_leaf(pvst, pvst_v_idx) || pvst.is_leaf(pvst_v_idx);

	return false;
}

void eval_vertex(const bd::VG &g, const pvst::Tree &pvst, pt::u32 pvst_v_idx,
		 std::vector<RoV> &rs)
{
	const pvst::VertexBase *pvst_v_ptr =
		pvst.get_vertex_const_ptr(pvst_v_idx);

	if (should_call(pvst, pvst_v_ptr, pvst_v_idx)) {
		RoV r{pvst_v_ptr};

		// get the set of walks for the RoV
		povu::genomics::graph::find_walks(g, r);

		// no walks found, skip this RoV
		if (r.get_walks().size() == 0)
			return;

		find_hidden(r);

		rs.push_back(std::move(r));
	}
}

/**
 * find walks in the graph based on the leaves of the pvst
 * initialize RoVs from flubbles
 */
std::vector<RoV> gen_rov(const std::vector<pvst::Tree> &pvsts, const bd::VG &g,
			 const core::config &app_config)
{

	// the set of RoVs to return
	std::vector<RoV> rs;
	rs.reserve(pvsts.size());

	// reuse msg buffer for progress bar
	DynamicProgress<ProgressBar> bars;
	std::string prog_msg;
	prog_msg.reserve(128);

	auto update_prog = [&](pt::idx_t i, pt::idx_t v_idx, pt::idx_t total,
			       std::size_t bar_idx)
	{
		prog_msg.clear();
		fmt::format_to(std::back_inserter(prog_msg),
			       "Generating RoVs for PVST {} ({}/{})", i + 1,
			       v_idx + 1, total);
		bars[bar_idx].set_option(option::PostfixText{prog_msg});
		bars[bar_idx].set_progress(v_idx + 1);
	};

	for (pt::idx_t i{}; i < pvsts.size(); i++) { // for each pvst
		const pvst::Tree &pvst = pvsts[i];
		// loop through each tree
		const pt::idx_t total = pvst.vtx_count();

		// reset progress bar
		ProgressBar bar;
		set_progress_bar_common_opts(&bar);
		std::size_t bar_idx = bars.push_back(bar);
		set_progress_bar_common_opts(&bar, pvst.vtx_count());

		for (pt::u32 v_idx{}; v_idx < pvst.vtx_count(); v_idx++) {
			// update progress bar
			if (app_config.show_progress())
				update_prog(i, v_idx, total, bar_idx);

			eval_vertex(g, pvst, v_idx, rs);

			if (app_config.show_progress())
				bars[bar_idx].mark_as_completed();
		}
	}

	return rs;
}

} // namespace povu::genomics::rov

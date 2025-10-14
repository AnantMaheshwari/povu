#ifndef POVU_GENOMICS_ROV_HPP
#define POVU_GENOMICS_ROV_HPP

#include <vector> // for vector

#include "povu/common/app.hpp"	     // for config
#include "povu/graph/bidirected.hpp" // for VG, bd
#include "povu/graph/pvst.hpp"	     // for Tree, VertexBase
#include "povu/graph/types.hpp"	     // for or_e, id_or_t, walk_t

namespace povu::genomics::rov
{
inline constexpr std::string_view MODULE = "povu::genomics::rov";
namespace pvst = povu::pvst;
namespace pgt = povu::types::graph;

enum class var_type_e {
	del, // deletion
	ins, // insertion
	sub, // substitution
	und  // undetermined
};

constexpr std::string_view to_string_view(var_type_e vt) noexcept
{
	switch (vt) {
	case var_type_e::del:
		return "DEL";
	case var_type_e::ins:
		return "INS";
	case var_type_e::sub:
		return "SUB";
	case var_type_e::und:
		return "UND";
	}
	// optional: handle out-of-range
	return "??";
}

inline std::ostream &operator<<(std::ostream &os, var_type_e vt)
{
	return os << to_string_view(vt);
}

struct rov_slice {
	pt::op_t<pt::u32> walk_idxs;
	// there's a one to one correspondence between slices and var_types
	std::vector<pt::slice_t> slices;
	std::vector<var_type_e> var_types;
};

/**
 * a collection of walks within a region of variation from start to end
 */
class RoV
{
	std::vector<pgt::walk_t> walks_;
	std::vector<rov_slice> extra_;
	const pvst::VertexBase *pvst_vtx;

public:
	// print when RoV is moved
	RoV(RoV &&other) noexcept
	    : walks_(std::move(other.walks_)), pvst_vtx(other.pvst_vtx)
	{
		other.pvst_vtx = nullptr;
	}

	RoV(const RoV &other) = delete;		   // disable copy constructor
	RoV &operator=(const RoV &other) = delete; // disable copy assignment
	RoV &
	operator=(RoV &&other) noexcept = delete; // disable move assignment

	// --------------
	// constructor(s)
	// --------------
	RoV(const pvst::VertexBase *v) : walks_(), pvst_vtx(v)
	{}

	// ---------
	// getter(s)
	// ---------

	[[nodiscard]]
	pt::idx_t walk_count() const
	{
		return this->walks_.size();
	}

	[[nodiscard]]
	const pvst::VertexBase *get_pvst_vtx() const
	{
		return this->pvst_vtx;
	}

	[[nodiscard]]
	const std::vector<pgt::walk_t> &get_walks() const
	{
		return this->walks_;
	}

	[[nodiscard]]
	std::vector<pgt::walk_t> &get_walks_mut()
	{
		return this->walks_;
	}

	[[nodiscard]]
	const std::vector<rov_slice> &get_extra() const
	{
		return this->extra_;
	}

	[[nodiscard]]
	const pgt::walk_t &get_walk(pt::idx_t idx) const
	{
		return this->walks_.at(idx);
	}

	// ---------
	// setter(s)
	// ---------

	void set_walks(std::vector<pgt::walk_t> &&walks)
	{
		this->walks_ = std::move(walks);
	}

	void add_extra(rov_slice &&rs)
	{
		this->extra_.emplace_back(std::move(rs));
	}

	// --------
	// other(s)
	// --------

	[[nodiscard]]
	std::string as_str() const
	{
		return this->pvst_vtx->as_str();
	}
};

/**
 * find walks in the graph based on the leaves of the pvst
 * initialize RoVs from flubbles
 */
std::vector<RoV> gen_rov(const std::vector<pvst::Tree> &pvsts, const bd::VG &g,
			 const core::config &app_config);

#ifdef TESTING
void find_hidden(RoV &r);
#endif

} // namespace povu::genomics::rov

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pgr = povu::genomics::rov;

#endif

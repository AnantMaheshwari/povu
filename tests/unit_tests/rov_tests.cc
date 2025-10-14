#include <gtest/gtest.h>
#include <vector>

#include "povu/genomics/graph.hpp"
#include "povu/genomics/rov.hpp"
#include "povu/graph/types.hpp"

namespace povu::unit_tests_rov
{
namespace pgg = povu::genomics::graph;

constexpr pgr::var_type_e ins = pgr::var_type_e::ins;
constexpr pgr::var_type_e del = pgr::var_type_e::del;
constexpr pgr::var_type_e sub = pgr::var_type_e::sub;

constexpr ptg::or_e fwd = ptg::or_e::forward;
constexpr ptg::or_e rev = ptg::or_e::reverse;

TEST(FindHidden, NarrowDown)
{
	pgr::RoV r{nullptr};
	ptg::walk_t w1{{1, fwd}, {2, fwd}, {3, fwd}, {5, fwd}};
	ptg::walk_t w2{{1, fwd}, {2, fwd}, {4, fwd}, {5, fwd}};
	r.set_walks({w1, w2});

	povu::genomics::rov::find_hidden(r);

	std::vector<pgr::rov_slice> expect_extra{
		{pt::op_t<pt::u32>{0, 1}, {{2, 1}}, {sub}},
		{pt::op_t<pt::u32>{1, 0}, {{2, 1}}, {sub}}};

	const auto &extra = r.get_extra();
	ASSERT_EQ(extra.size(), expect_extra.size());
	for (pt::idx_t i{}; i < extra.size(); i++) {
		EXPECT_EQ(extra[i].walk_idxs, expect_extra[i].walk_idxs);
		EXPECT_EQ(extra[i].slices, expect_extra[i].slices);
		EXPECT_EQ(extra[i].var_types, expect_extra[i].var_types);
	}
}

TEST(FindHidden, IndelsAndNarrowDown)
{
	pgr::RoV r{nullptr};
	ptg::walk_t w1{{1, fwd}, {2, fwd}, {3, fwd}, {5, fwd}};
	ptg::walk_t w2{{1, fwd}, {3, fwd}, {5, fwd}};
	r.set_walks({w1, w2});

	povu::genomics::rov::find_hidden(r);

	std::vector<pgr::rov_slice> expect_extra{
		{pt::op_t<pt::u32>{0, 1}, {{1, 1}}, {ins}},
		{pt::op_t<pt::u32>{1, 0}, {{0, 1}}, {del}}};

	const auto &extra = r.get_extra();
	ASSERT_EQ(extra.size(), expect_extra.size());
	for (pt::idx_t i{}; i < extra.size(); i++) {
		EXPECT_EQ(extra[i].walk_idxs, expect_extra[i].walk_idxs);
		EXPECT_EQ(extra[i].slices, expect_extra[i].slices);
		EXPECT_EQ(extra[i].var_types, expect_extra[i].var_types);
	}
}

TEST(FindHidden, IndelsAndMismatchNarrowDown)
{
	pgr::RoV r{nullptr};
	ptg::walk_t w1{{1, fwd}, {2, fwd}, {4, fwd}, {5, fwd}};
	ptg::walk_t w2{{1, fwd}, {3, fwd}, {5, fwd}};
	r.set_walks({w1, w2});

	povu::genomics::rov::find_hidden(r);

	std::vector<pgr::rov_slice> expect_extra{
		{pt::op_t<pt::u32>{0, 1}, {{1, 2}}, {sub}},
		{pt::op_t<pt::u32>{1, 0}, {{1, 1}}, {sub}}};

	const auto &extra = r.get_extra();
	ASSERT_EQ(extra.size(), expect_extra.size());
	for (pt::idx_t i{}; i < extra.size(); i++) {
		EXPECT_EQ(extra[i].walk_idxs, expect_extra[i].walk_idxs);
		EXPECT_EQ(extra[i].slices, expect_extra[i].slices);
		EXPECT_EQ(extra[i].var_types, expect_extra[i].var_types);
	}
}

} // namespace povu::unit_tests_rov

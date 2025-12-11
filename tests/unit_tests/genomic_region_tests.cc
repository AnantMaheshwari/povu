#include <gtest/gtest.h>
#include <optional>
#include <string>

#include "povu/variation/rov.hpp"

namespace pvr = povu::var::rov;

// Test parse_genomic_region function
TEST(GenomicRegionTest, ParseValidRegion)
{
	auto region = pvr::parse_genomic_region("chr17:43044294-43125482");
	ASSERT_TRUE(region.has_value());
	EXPECT_EQ(region->ref_name, "chr17");
	EXPECT_EQ(region->start, 43044294);
	EXPECT_EQ(region->end, 43125482);
}

TEST(GenomicRegionTest, ParseValidRegionWithPrefix)
{
	auto region = pvr::parse_genomic_region("GRCh38#chr17:100-200");
	ASSERT_TRUE(region.has_value());
	EXPECT_EQ(region->ref_name, "GRCh38#chr17");
	EXPECT_EQ(region->start, 100);
	EXPECT_EQ(region->end, 200);
}

TEST(GenomicRegionTest, ParseInvalidFormatMissingColon)
{
	auto region = pvr::parse_genomic_region("chr17_100-200");
	EXPECT_FALSE(region.has_value());
}

TEST(GenomicRegionTest, ParseInvalidFormatMissingDash)
{
	auto region = pvr::parse_genomic_region("chr17:100");
	EXPECT_FALSE(region.has_value());
}

TEST(GenomicRegionTest, ParseInvalidFormatNonNumericStart)
{
	auto region = pvr::parse_genomic_region("chr17:abc-200");
	EXPECT_FALSE(region.has_value());
}

TEST(GenomicRegionTest, ParseInvalidFormatNonNumericEnd)
{
	auto region = pvr::parse_genomic_region("chr17:100-xyz");
	EXPECT_FALSE(region.has_value());
}

TEST(GenomicRegionTest, ParseInvalidFormatEmptyStart)
{
	auto region = pvr::parse_genomic_region("chr17:-200");
	EXPECT_FALSE(region.has_value());
}

TEST(GenomicRegionTest, ParseInvalidFormatEmptyEnd)
{
	auto region = pvr::parse_genomic_region("chr17:100-");
	EXPECT_FALSE(region.has_value());
}

TEST(GenomicRegionTest, ParseInvalidFormatStartGreaterThanEnd)
{
	auto region = pvr::parse_genomic_region("chr17:200-100");
	EXPECT_FALSE(region.has_value());
}

TEST(GenomicRegionTest, ParseInvalidFormatStartEqualsEnd)
{
	auto region = pvr::parse_genomic_region("chr17:100-100");
	EXPECT_FALSE(region.has_value());
}

TEST(GenomicRegionTest, ParseInvalidFormatEmptyRefName)
{
	auto region = pvr::parse_genomic_region(":100-200");
	EXPECT_FALSE(region.has_value());
}

TEST(GenomicRegionTest, GenomicRegionIsValid)
{
	pvr::genomic_region valid_region("chr1", 100, 200);
	EXPECT_TRUE(valid_region.is_valid());

	pvr::genomic_region invalid_region1("", 100, 200);
	EXPECT_FALSE(invalid_region1.is_valid());

	pvr::genomic_region invalid_region2("chr1", 200, 100);
	EXPECT_FALSE(invalid_region2.is_valid());

	pvr::genomic_region invalid_region3("chr1", 100, 100);
	EXPECT_FALSE(invalid_region3.is_valid());
}

TEST(GenomicRegionTest, ParseLargeNumbers)
{
	auto region = pvr::parse_genomic_region("chr1:1000000000-2000000000");
	ASSERT_TRUE(region.has_value());
	EXPECT_EQ(region->ref_name, "chr1");
	EXPECT_EQ(region->start, 1000000000);
	EXPECT_EQ(region->end, 2000000000);
}

TEST(GenomicRegionTest, ParseZeroStart)
{
	auto region = pvr::parse_genomic_region("chr1:0-100");
	ASSERT_TRUE(region.has_value());
	EXPECT_EQ(region->ref_name, "chr1");
	EXPECT_EQ(region->start, 0);
	EXPECT_EQ(region->end, 100);
}

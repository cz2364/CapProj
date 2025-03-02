#include "defs.h"
#include<string>
#include<vector>
#include <gtest/gtest.h>

TEST(MyLibTest, AddTest) {
    std::vector<std::string> vn = {"x", "y"};
    std::vector<long> vi = {0, 1};
    int cdim  = 2;
    std::string cn = "test_cone";
    StdSoc test_cone = StdSoc(cn, cdim, vn, vi);
    EXPECT_EQ(test_cone.get_name(), cn);
    EXPECT_EQ(test_cone.get_dim(), cdim);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


#include "defs.h"
#include "utils.h"
#include<string>
#include<vector>
#include <gtest/gtest.h>

TEST(DefsTest, StdScoCreationTest) {
    std::vector<std::string> vn = {"x", "y"};
    std::vector<long> vi = {0, 1};
    int cdim  = 2;
    std::string cn = "test_cone";
    StdSoc test_cone = StdSoc(cn, cdim, vn, vi);
    EXPECT_EQ(test_cone.get_name(), cn);
    EXPECT_EQ(test_cone.get_dim(), cdim);
}

TEST(UtilsTest, ClipTest){
    std::vector<double> ub = {3.0, 2.0, 5.0};
    std::vector<double> lb = {-2.1, 1.5, -5.0};
    std::vector<double> x1 = {1.5, 2.2, 5.1};
    std::vector<double> px1 = {1.5, 2.0, 5.0};
    std::vector<double> x2 = {1.4, 2.1, 0.0};
    std::vector<double> px2 = {1.5, 2.0, 0.0};
    std::vector<double> x3 = {1.4, 2.1, 0.0};

    //std::vector<double> tx1 = box_proj(x1, ub, lb);
    //for(int i = 0; i < tx1.size(); i++){
    //    EXPECT_EQ(tx1[i], px1[i]);
    //}

    // std::vector<double> tx2 = box_proj(x2, ub, lb);
    //for(int i = 0; i < tx2.size(); i++){
    //    EXPECT_EQ(tx2[i], px2[i]);
    // }


}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


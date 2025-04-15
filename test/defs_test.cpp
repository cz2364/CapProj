#include "defs.h"
#include "utils.h"
#include<string>
#include<vector>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <gtest/gtest.h>
#include "absl/log/log.h"
#include "absl/log/check.h"


TEST(UtilsTest, ClipTest){
    std::vector<double> ub = {3.0, 2.0, 5.0};
    std::vector<double> lb = {-2.1, 1.5, -5.0};
    std::vector<double> x1 = {1.5, 2.2, 5.1};
    std::vector<double> px1 = {1.5, 2.0, 5.0};
    std::vector<double> x2 = {1.4, 2.1, 0.0};
    std::vector<double> px2 = {1.4, 2.0, 0.0};
    std::vector<double> x3 = {1.4, 2.1, 0.0};

    std::vector<double> tx1 = UTILS_H::box_proj(x1, ub, lb);
    for(int i = 0; i < tx1.size(); i++){
        EXPECT_EQ(tx1[i], px1[i]);
    }

    std::vector<double> tx2 = UTILS_H::box_proj(x2, ub, lb);
    for(int i = 0; i < tx2.size(); i++){
        EXPECT_EQ(tx2[i], px2[i]);
    }

}

TEST(UtilsTest, VarNamesTest){
    std::string s1 = "x";
    std::string s2 = "x__COPY_1";
    std::string s3 = "x__COPY_";
    std::string s4 = "x__COPY_a";
    std::string s5 = "x__COPY_11";
    EXPECT_EQ(generate_variable_name(s1), "x__COPY_1");
    EXPECT_EQ(generate_variable_name(s2), "x__COPY_2");
    EXPECT_EQ(generate_variable_name(s3), "x__COPY___COPY_1");
    EXPECT_EQ(generate_variable_name(s4), "x__COPY_a__COPY_1");
    EXPECT_EQ(generate_variable_name(s5), "x__COPY_12");
}


TEST(DefsTest, EigenMatrixTest){
    Eigen::MatrixXf m(2, 3);
    m << 1.0, 2.0, 3.0,
         4.0, 5.0, 6.0;
    EXPECT_EQ(m.rows(), 2);
    EXPECT_EQ(m.cols(), 3);

}




int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    LOG(INFO) << "Testing";

    return RUN_ALL_TESTS();
}


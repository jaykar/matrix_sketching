#include <limits.h>
#include "../../armadillo/sk_arm.hpp"
#include "gtest/gtest.h"
#include <iostream>

using namespace std; 
TEST(construct_arm_interface, TRIVIAL){
    for(int i=0; i<1000; i++){
        auto A = sk_arm();
    }
    EXPECT_EQ(1,1);
}


TEST(construct_arm_interface_num, TRIVIAL){
    for(int i=0; i<1000; i++){
        auto A = sk_arm(20, 20);
    }
    EXPECT_EQ(1,1);
}

TEST(construct_arm_interface_mat_copy, TRIVIAL){
    auto B = sk_arm(30,30); 
    for(int i=0; i<1000; i++){
        auto A = B;
    }
    EXPECT_EQ(1,1);
}


TEST(construct_arm_interface_mat_move, TRIVIAL){
    auto B = sk_arm(); 
    B.rand_n(30,30); 
    auto C = B; 
    auto A = std::move(B); 
    for(int i=0; i<1000; i++){
        B = C; 
        A = std::move(B); 
    }
    EXPECT_EQ(1,1);
}

int main(int argc, char **argv){
    cout << "this is to test for memory leaks for constructors" << endl; 
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

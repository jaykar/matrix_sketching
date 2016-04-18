#include <limits.h>
#include "../../Arm_interface.hpp"
#include "gtest/gtest.h"
#include <iostream>

TEST(construct_arm_interface, TRIVIAL){
    std::cout << "this is to test for memory leaks when passing default" << std::endl;  
    for(int i=0; i<1000; i++){
        auto A = Armadillo_Matrix(); 
    }

    //ignore this it's just to regist test
    EXPECT_EQ(1,1); 
}

TEST(construct_arm_interface_row_col, TRIVIAL){
    std::cout << "this is to test for memory leaks when passing row, col" << std::endl;  
    for(int i=0; i<1000; i++){
        auto A = Armadillo_Matrix(30, 30); 
    }

    //ignore this it's just to register test
    EXPECT_EQ(1,1); 
}

Armadillo_Matrix f(){
    auto A = Armadillo_Matrix(); 
    cout << A.data() << endl; 
    return A.rand_n(3,4, 0, 1.0); 
}
int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS(); 
}

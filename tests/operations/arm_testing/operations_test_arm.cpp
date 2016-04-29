#include <iostream>
#include "operations.hpp"
//#include "boost.hpp"
#include "armadillo.hpp"
#include <string>
using namespace std;

//namespace bnu = boost::numeric::ublas;

void ksvd_test(){
    string x_label = "data_image";
    cout << "loading x train" << endl;
    sketchy::armadillo x_train = sketchy::armadillo(x_label);

    sketchy::armadillo U;
    sketchy::armadillo S;
    sketchy::armadillo V;
    sketchy::ops::k_svd(x_train, U, S, V, 5, 20);

    U.save("u");
    S.save("s");
    V.save("v");
}
    /*
    cout << "loading y train" << endl;
    string y_label = "../../../../mnist/regularized_label_small";
    sketchy::armadillo y_train = sketchy::armadillo(y_label);

    cout << "solving for b" << endl;
    sketchy::armadillo b = x_train.solve_x(y_train);

    cout << "solving for xtilde" << endl;
    sketchy::armadillo x_tilde = sketchy::ops::lin_regress<sketchy::armadillo>(x_train, y_train, 600);

    //sketchy::armadillo y_pred = x_train.mult(x_tilde);
    //sketchy::armadillo y_pred2 = x_train.mult(b);

    b.save("b");
    x_tilde.save("x_tilde");
    */
    //cout << (y_pred.subtract(y_pred2)).accumulate()  << endl;

    return 0;
}

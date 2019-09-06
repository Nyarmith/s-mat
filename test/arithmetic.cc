#include "s-mat.hh"
#include "util.hh"
#include <cassert>
#include <iostream>

using namespace SMat;

int main()
{
    Vec<int,3> A{1,2,3};
    Vec<int,3> B{4,5,6};
    assert(A*B == 32);

    Vec<int,3> AB{5,7,9};
    assert( A+B == AB );

    Matrix<int,4,3> C {1, 5, 9,
                       2, 6, 10,
                       3, 7, 11,
                       4, 8, 12};

    Matrix<int,3,4> D {1,  2,  3,  4,
                       5,  6,  7,  8,
                       9,  10, 11, 12};

	Matrix<int,4,4> E{107, 122, 137, 152,
                      122, 140, 158, 176,
                      137, 158, 179, 200,
                      152, 176, 200, 224};

	assert(C*D == E);

    Matrix<double,4,3> F {.5,   .25,    .125,
                          .0652, .03125, .015625,
                          .3,    .09,    .027,
                          .0081, .00243, .000729};

   Matrix<double,3,4> F_prime  {0.5,   0.0652,   0.3,   0.0081,
                                0.25,  0.03125,  0.09,  0.00243,
                                0.125, 0.015625, 0.027, 0.000729};
   


   Matrix<double,4,4> F_prod {0.328125,    0.042365625, 0.175875, 0.004748625,
                              0.042365625, 0.005471743, 0.022794375, 0.000615448,
                              0.175875,    0.022794375, 0.098829,    0.002668383,
                              0.004748625, 0.000615448, 0.002668383, 0.000072046};

   assert(is_close(F*F_prime,F_prod,0.0000001));

   Matrix<int,3,3> G {
       1, 5, 9,
       2, 6, 10,
       3, 7, 11};

   Matrix<int,3,3> G_3 {
       708, 1884, 3060,
       840, 2232, 3624,
       972, 2580, 4188};

   auto cG_3 = pow(G,3);
   assert(cG_3 == G_3);
   
   
   assert( (G[0] == Vec<int,3>{1,5,9}) );
   assert(G[0][2] == 9);
   assert(G[2][0] == 3);


   Matrix<int,3,3> H {1, 1, 1,
                      1, 1, 1,
                      1, 1, 1}; //equivalent to explcitl c-tor w/ arg 1

   Matrix<int,3,3> G_m_H {
       0,4,8,
       1,5,9,
       2,6,10};

   assert(G-H == G_m_H);
   assert(G_m_H + H == G);
   assert(G != H);

   Matrix<int,3,3> Threes(3);
   assert(H*3 == Threes);
   assert(3*H == Threes);
   assert(H == Threes/3);
}


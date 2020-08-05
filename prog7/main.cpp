//
//  main.cpp
//  yuw14_prog7
//
//  Created by 王畬 on 2020/7/18.
//  Copyright © 2020 Yu Wang. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include "/Users/wangyu/Downloads/IE531/prog7/newmat11/newmatio.h"
#include "/Users/wangyu/Downloads/IE531/prog7/newmat11/newmat.h"
#include "/Users/wangyu/Downloads/IE531/prog7/newmat11/newmatap.h"
#include "/Users/wangyu/Downloads/IE531/prog7/newmat11/include.h"

using namespace std;

Matrix Sketching(double epsilon, Matrix input){
    int rows = input.nrows();
    int Bcols = ceil(2.0/epsilon);
    Matrix U,V;
    DiagonalMatrix Xigma;
    if(Bcols>=rows){
        SVD(input*input.t(), Xigma, U, V);
        
        Matrix B(rows,rows);
        B = 0.0;
        Matrix Xigma_(rows,rows);
        Xigma_ = 0.0;
        
        for(int i=1;i<=rows;++i){
            Xigma_(i,i)=sqrt(Xigma(i,i));
        }
        B=U*Xigma_*V.t();
        
        return B;
    }
    else{
        Matrix B(rows,Bcols);
        B = 0.0;
        Matrix Xigma_(Bcols,Bcols);
        Xigma_ = 0.0;
        
        ColumnVector Empty_Column(rows);
        Empty_Column = 0.0;
        
        for(int i=1;i<=input.ncols();++i){
            for(int j=1;j<=B.ncols();++j){
                if(B.Column(j)!=Empty_Column){
                    SVD(B, Xigma, U, V);
                    
                    IdentityMatrix I_l(Xigma.nrows());
                    double sigma_i = Xigma(Bcols)*Xigma(Bcols);// the least singular value of B
                    Xigma = Xigma*Xigma.t()-sigma_i*I_l;
                    for(int k=1;k<=Bcols;++k){
                        Xigma_(k,k)=sqrt(Xigma(k,k));
                    }
                    B = U*Xigma_;
                }
                else{
                    B.Column(j)=input.Column(i);
                    break;
                }
            }
        }
        return B;
    }
}

double frobeniusNorm(Matrix input){
    Matrix m = input * input.t();
    int r = m.nrows();
    int c = m.ncols();
    double f_n = 0.0;
    for(int i=1;i<=r;i++){
        for(int j=1;j<=c;j++){
            f_n = f_n + m(i,j) * m(i,j);
        }
    }
    return sqrt(f_n);
}

int main(int argc, const char * argv[]) {
    double epsilon;
    int rows,columns;
    sscanf(argv[1], "%d", &rows);
    sscanf(argv[2], "%d", &columns);
    sscanf(argv[3], "%lf", &epsilon);
    ifstream inputFile(argv[4]);
    ofstream outputFile(argv[5]);
    
    Matrix input_matrix(rows,columns);
    for(int i=1;i<=rows;i++){
        for(int j=1;j<=columns;j++){
            inputFile>>input_matrix(i,j);
        }
    }
    Matrix output_matrix(rows,min((int)ceil(2/epsilon),rows));
    output_matrix = Sketching(epsilon, input_matrix);
    
    double original_FN = frobeniusNorm(input_matrix);
    double new_FN = frobeniusNorm(output_matrix);
    
    cout<<"Edo Liberty's Matrix Sketching Algorithm"<<endl;
    cout<<"----------------------------------------"<<endl;
    cout<<"Original Data-Matrix has "<<rows<<"-rows & "<<columns<< "-cols"<<endl;
    cout<<"Epsilon = "<<epsilon<<" (i.e. max. of "<<100*epsilon<<"% reduction of Frobenius-Norm of the Sketch Matrix)"<< endl;
    cout<<"Input File = "<<argv[4]<<endl;
    cout<<"Frobenius Norm of the ("<<input_matrix.nrows()<<" x "<< input_matrix.ncols()<<") Data Matrix = "<<original_FN<<endl;
    cout<<"Frobenius Norm of the ("<<output_matrix.nrows()<<" x "<<output_matrix.ncols()<<") Sketch Matrix = "<<new_FN<<endl;
    cout<<"Change in Frobenius-Norm between Sketch & Original  = "<<setprecision(3)<<(new_FN-original_FN)*100/original_FN<<"%"<<endl;
    cout<<"File '"<<argv[5]<<"' contains a ("<<output_matrix.nrows()<<" x "<<output_matrix.ncols()<<") Matrix-Sketch"<<endl;

    outputFile<<output_matrix;
    
}

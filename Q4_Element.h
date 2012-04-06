#ifndef Q4_ELEMENT_H
#define Q4_ELEMENT_H
#include <vector>
#include <string>
#include "Node.h"
#include "ElasticMaterial.h"
#include <cmath>
#include <Eigen/Dense>
#include <iostream>

using namespace std;
using namespace Eigen;

class Q4_Element
{
    public:

        Q4_Element();
        Q4_Element(int,double,double,vector <Node>);
        void SetId(int);
        int GetId(void);
        void Setth(double);
        double Getth(void);
        void SetGama(double);
        double GetGama(void);
        void SetMat(ElasticMaterial);
        ElasticMaterial GetMat(void);
        void SetNodalObj(vector <Node>);
        double Calc_ShapeFunction(int,int,double,double);
        double Calc_DiffN(int,int,double,double,string);
        MatrixXd Calc_BMatrix(double,double);
        void SetDMatrix(string);
        void Calc_LSM();
        MatrixXd Get_LSM();
        void SetU(MatrixXd);
        MatrixXd Get_Ep();
        MatrixXd Get_Sigma();
        MatrixXd Get_PSigma(string);
        MatrixXd Get_PStrain(string);

        virtual ~Q4_Element();

    protected:

    private:

        int id;
        vector <Node> NodeObj;
        ElasticMaterial mat;
        double E,v,gama,th;
        double gpt [2]={-1/sqrt(3),1/sqrt(3)};
        double wgpt[2]={1,1};
        MatrixXd LSM=MatrixXd::Zero(8,8);
        MatrixXd B=MatrixXd::Zero(3,8);
        MatrixXd D=MatrixXd::Zero(3,3);
        MatrixXd U=MatrixXd::Zero(8,1);
        MatrixXd Ep=MatrixXd::Zero(3,4);
        MatrixXd Sigma=MatrixXd::Zero(3,4);
        MatrixXd PSigma=MatrixXd::Zero(3,4);
        MatrixXd PStrain=MatrixXd::Zero(3,4);
};

#endif // Q4_ELEMENT_H

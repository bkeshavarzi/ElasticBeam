#ifndef Q9_ELEMENT_H
#define Q9_ELEMENT_H
#include <iostream>
#include <Eigen/Dense>
#include "ElasticMaterial.h"
#include <vector>
#include <string>
#include "Node.h"
#include <cmath>

using namespace std;
using namespace Eigen;

class Q9_Element
{
    public:
        Q9_Element();
        Q9_Element(int,double,vector <Node>,ElasticMaterial);
        void SetId(int);
        int GetId(void);
        void Setth(double);
        double Getth(void);
        void SetGama(double);
        double GetGama(void);
        void SetMat(ElasticMaterial);
        ElasticMaterial GetMat(void);
        void SetNodalObj(vector <Node>);
        vector <Node> GetNodalObj(void);
        void SetElemParam(double []);
        double f_kesi(int,double);
        double f_eta(int,double);
        double Calc_ShapeFunction(int,int,double,double);
        double Calc_DiffN(int,int,double,double,string);
        MatrixXd Calc_BMatrix(double,double);
        void SetDMatrix(string);
        void Setlocalcord(void);
        MatrixXd Getlocalcord(void);
        MatrixXd CalcJacobian(double,double);
        MatrixXd CalcInvJacobian(double,double);
        double CalcDetJacobian(double,double);
        MatrixXd GetLocalCord();
        MatrixXd Calc_LSM();
        void SetU(MatrixXd);
        MatrixXd GetU();
        void Set_Ep();
        void Set_Sigma();
        void Set_PSigma(string);
        void Set_PStrain(string);
        MatrixXd Get_Ep();
        MatrixXd Get_Sigma();
        MatrixXd Get_PSigma();
        MatrixXd Get_PStrain();
        friend int SignFunction(double);
        virtual ~Q9_Element();

    protected:

    private:

        int id;
        vector <Node> NodeObj;
        ElasticMaterial mat;
        double E,v,gama,th,a,b,xc,yc;
        double gpt [3]={-1*sqrt(0.6),0,sqrt(0.6)};
        double wgpt[3]={5/9,8/9,5/9};
        MatrixXd LCord=MatrixXd::Zero(2,9);
        MatrixXd LSM=MatrixXd::Zero(18,18);
        MatrixXd D=MatrixXd::Zero(3,3);
        MatrixXd U=MatrixXd::Zero(18,1);
        MatrixXd Ep=MatrixXd::Zero(3,9);
        MatrixXd Sigma=MatrixXd::Zero(3,9);
        MatrixXd PSigma=MatrixXd::Zero(3,9);
        MatrixXd PStrain=MatrixXd::Zero(3,9);
};

#endif // Q9_ELEMENT_H

#ifndef Q8_ELEMENT_H
#define Q8_ELEMENT_H

#include "Q4_Element.h"
#include "Q9_Element.h"
//WATCH OUT FOR NODE NUMBERING IN MESH DATA AND CONSIDERED NODE NUMERING IN ELEMENT CLASS

class Q8_Element : public Q4_Element,Q9_Element
{
    public:

        Q8_Element();
        Q8_Element(int,double,vector <Node>);
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
        friend int SignFunction(double);
        virtual ~Q8_Element();

    protected:

    private:

        int id;
        vector <Node> NodeObj;
        ElasticMaterial mat;
        double E,v,gama,th;
        double gpt [3]={-1*sqrt(0.6),0,sqrt(0.6)};
        double wgpt[3]={5/9,8/9,5/9};
        MatrixXd LSM=MatrixXd::Zero(16,16);
        MatrixXd B=MatrixXd::Zero(3,18);
        MatrixXd D=MatrixXd::Zero(3,3);
        MatrixXd U=MatrixXd::Zero(16,1);
        MatrixXd Ep=MatrixXd::Zero(3,8);
        MatrixXd Sigma=MatrixXd::Zero(3,8);
        MatrixXd PSigma=MatrixXd::Zero(3,8);
        MatrixXd PStrain=MatrixXd::Zero(3,8);
};

#endif // Q8_ELEMENT_H

#ifndef CSE_ELEMENT_H
#define CSE_ELEMENT_H
#include "Node.h"
#include <ElasticMaterial.h>
#include <Eigen/Cholesky>
#include <cmath>
#include <string>

class CSE_Element
{
    public:
        CSE_Element();
        CSE_Element(int,double,vector <Node>,ElasticMaterial); //id,thickness,nodal objects,material
        void SetId(int);
        int GetId(void);
        void Setth(double);
        double Getth(void);
        void SetNodeObj(Node,Node,Node);
        vector <Node> GetNodeObj(void);
        void SetDof();
        void SetCMatrix();
        MatrixXd GetCMatrix();
        void SetA();
        double GetA();
        void SetBMatrix();
        MatrixXd GetBMatrix(void);
        void SetDMatrix_PlaneStress();
        void SetDMatrix_PlaneStrain();
        MatrixXd GetDMatrix();
        void SetLSM();
        MatrixXd GetLSM();
        void SetU(MatrixXd);
        void CalculateStrainTensor();
        MatrixXd GetStrainTensor();
        void CalculateStressTensor();
        MatrixXd GetStressTensor();
        void CalculatePrinStress(string);
        MatrixXd GetPrinStress();
        void CalculatePrinStrain();
        MatrixXd GetPrinStrain();
        virtual ~CSE_Element();

    protected:

    private:
        int id;
        vector <Node> NodeObj;
        int * ldof;
        ElasticMaterial mat;
        double E,v,gama,th,A;
        MatrixXd C=MatrixXd::Zero(3,3);
        MatrixXd LSM=MatrixXd::Zero(6,6);
        MatrixXd B=MatrixXd::Zero(3,6);
        MatrixXd D=MatrixXd::Zero(3,3);
        MatrixXd U=MatrixXd::Zero(6,1);
        MatrixXd Ep=MatrixXd::Zero(3,1);
        MatrixXd Sigma=MatrixXd::Zero(3,1);
        MatrixXd PSigma=MatrixXd::Zero(3,1);
        MatrixXd PStrain=MatrixXd::Zero(3,1);
};

#endif // CSE_ELEMENT_H

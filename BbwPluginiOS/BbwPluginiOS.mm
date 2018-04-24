//
//  BbwPluginiOS.mm
//  BbwPluginiOS
//
//  Created by Zeratul Phillip on 18/04/2018.
//  Copyright © 2018 Zeratul Phillip. All rights reserved.
//

#import "BbwPluginiOS.h"
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMatrixXd;

extern "C" {
    
    int BbwPluginIOS(float* ptr_V, int vCount, int oriVCount,
                     int* ptr_F, int fCount,
                     float* ptr_C, int cCount,
                     int* ptr_BE, int beCount,
                     double* ptr_W)
    {
        std::cout << "Testing BbwPluginIOS" << std::endl;

        // V : Vertices : std::vector<std::vector<double> > : numOfVertices * double[3] -> Convert from float[][2] to double [][3]
        // F : Triangle : Index std::vector<std::vector<int> > : numOfTriangles * int[3] -> Use directly
        /// dont need T : Tetrahedra : Index std::vector<std::vector<int> > : numOfTetrahedra * int[4]
        // C : TGF Vertex : std::vector<std::vector<double> > : numOfVertices * double[3] -> Convert from float[][2] to double [][3]
        // BE : TGF Edges : std::vector<std::vector<int> > : numOfEdges * int[2] -> Use directly
        // W : Weights Output

        double array_V[vCount][3];
        for (int i = 0; i < vCount; i++)
        {
            array_V[i][0] = (double)ptr_V[i * 2];
            array_V[i][1] = (double)ptr_V[i * 2 + 1];
            array_V[i][2] = (double)0;
        }
        double array_C[cCount][3];
        for (int i = 0; i < vCount; i++)
        {
            array_C[i][0] = (double)ptr_C[i * 2];
            array_C[i][1] = (double)ptr_C[i * 2 + 1];
            array_C[i][2] = (double)0;
        }
        
    
        Eigen::MatrixXd V, C, W;
        Eigen::MatrixXi T, F, BE;
        // List of boundary indices (aka fixed value indices into VV)
        Eigen::VectorXi b;
        // List of boundary conditions of each weight function
        Eigen::MatrixXd bc;

        V = Eigen::Map<Eigen::MatrixXd>(&array_V[0][0], vCount, 3);
        F = Eigen::Map<Eigen::MatrixXi>(ptr_F, fCount, 3);
        C = Eigen::Map<Eigen::MatrixXd>(&array_C[0][0], cCount, 3);
        BE = Eigen::Map<Eigen::MatrixXi>(ptr_BE, beCount, 2);
        
        
        
        ///
//        std::string sep = "\n----------------------------------------\n";
//        Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
//        std::cout << V.format(CleanFmt) << sep;
//        std::cout << F.format(CleanFmt) << sep;
//        std::cout << C.format(CleanFmt) << sep;
//        std::cout << BE.format(CleanFmt) << sep;
        ///
        std::cout << "V[" << V.rows() << ", " << V.cols() << "] " << vCount << std::endl;
        std::cout << "F[" << F.rows() << ", " << F.cols() << "] " << fCount << std::endl;
        std::cout << "C[" << C.rows() << ", " << C.cols() << "] " << cCount << std::endl;
        std::cout << "BE[" << BE.rows() << ", " << BE.cols() << "] " << beCount << std::endl;
        
        
        igl::boundary_conditions(V,T,C,Eigen::VectorXi(),BE,Eigen::MatrixXi(),b,bc);

        
        std::cout << "Flag 1" << std::endl;
        
        
        // compute BBW weights matrix
        igl::BBWData bbw_data;
        // only a few iterations for sake of demo
        bbw_data.active_set_params.max_iter = 8;
        bbw_data.verbosity = 2;
        //if(!igl::bbw(V,T,b,bc,bbw_data,W))
        if(!igl::bbw(V,F,b,bc,bbw_data,W))
        {
            return EXIT_FAILURE;
        }
        
        std::cout << "Flag 2" << std::endl;
        
        // Normalize weights to sum to one
        igl::normalize_row_sums(W,W);
        
        Eigen::Map<Eigen::MatrixXd>(ptr_W, W.rows(), W.cols()) = W;
        
        return 0;
    }
    
//    int BbwPluginIOS(float* ptr_V, int vCount, int oriVCount,
//                     int* ptr_F, int fCount,
//                     float* ptr_C, int cCount,
//                     int* ptr_BE, int beCount,
//                     float* ptr_W)
//    {
//        std::cout << "Testing BbwPluginIOS" << std::endl;
//
//        for (int i = 0; i < vCount * 2; i++)
//            std::cout << "ptr_V[" << i << "] = " << ptr_V[i] << std::endl;
//        std::cout << "oriVCount = " << oriVCount << std::endl;
//        for (int i = 0; i < fCount; i++)
//            std::cout << "ptr_F[" << i << "] = " << ptr_F[i] << std::endl;
//        for (int i = 0; i < cCount * 2; i++)
//            std::cout << "ptr_C[" << i << "] = " << ptr_C[i] << std::endl;
//        for (int i = 0; i < beCount * 2; i++)
//            std::cout << "ptr_BE[" << i << "] = " << ptr_BE[i] << std::endl;
//
//        for (int i = 0; i < cCount * oriVCount; i++)
//                ptr_W[i] = i*10;
//
//        for (int i = 0; i < cCount * oriVCount; i++)
//            std::cout << "ptr_W[" << i << "] = " << ptr_W[i] << std::endl;
//
//        return 0;
//    }
}

//
//  BbwPluginiOS.mm
//  BbwPluginiOS
//
//  Created by Zeratul Phillip on 18/04/2018.
//  Copyright © 2018 Zeratul Phillip. All rights reserved.
//

#import "BbwPluginiOS.h"

extern "C" {
    
    int BbwPluginIOS(float* ptr_V, int vCount, int oriVCount,
                     int* ptr_F, int fCount,
                     float* ptr_C, int cCount,
                     int* ptr_BE, int beCount,
                     double* ptr_W)
    {
        std::cout << "BbwPluginIOS" << std::endl;
        
        // A little confused with colMajor and rowMajor, so maybe not correct here.
        // V  (Vertices)             : float[vCount, 2](rowMajor) -> double [vCount * 3](colMajor) -> Matrix<double, vCount, 3, colMajor>
        // F  (Triangle Indices)     : int[fCount * 3](rowMajor)  -> int[fCount * 3](colMajor)     -> Matrix<int, fCount, 3, colMajor>
        // C  (Bone Point Positions) : float[cCount, 2](rowMajor) -> double[cCount * 3](colMajor)  -> Matrix<double, cCount, 3, colMajor>
        // BE (Bone Edge Indices)    : int[beCount, 2](rowMajor)  -> int[beCount * 2](colMajor)    -> Matrix<int, beCount, 2, colMajor>
        // W  (Output Weights)       : double[beCount, oriVCount] <- double[beCount * oriVCount]   <- Matrix<double, oriCount, beCount, colMajor>
        
        double array_V[vCount * 3];
        for (int i = 0; i < vCount; i++)
        {
            array_V[i] = (double)(ptr_V[i * 2]);
            array_V[i + vCount] = (double)(ptr_V[i * 2 + 1]);
            array_V[i + vCount * 2] = (double)0;
        }
        int array_F[fCount * 3];
        for (int i = 0; i < fCount; i++)
        {
            array_F[i] = ptr_F[i * 3];
            array_F[i + fCount] = ptr_F[i * 3 + 1];
            array_F[i + fCount * 2] = ptr_F[i * 3 + 2];
        }
        double array_C[cCount * 3];
        for (int i = 0; i < cCount; i++)
        {
            array_C[i] = (double)(ptr_C[i * 2]);
            array_C[i + cCount] = (double)(ptr_C[i * 2 + 1]);
            array_C[i + cCount * 2] = (double)0;
        }
        int array_BE[beCount * 2];
        for (int i = 0; i < beCount; i++)
        {
            array_BE[i] = ptr_BE[i * 2];
            array_BE[i + beCount] = ptr_BE[i * 2 + 1];
        }
        
        Eigen::MatrixXd V, C, W;
        Eigen::MatrixXi T, F, BE;
        // List of boundary indices (aka fixed value indices into VV)
        Eigen::VectorXi b;
        // List of boundary conditions of each weight function
        Eigen::MatrixXd bc;

        V = Eigen::Map<Eigen::MatrixXd>(&array_V[0], vCount, 3);
        F = Eigen::Map<Eigen::MatrixXi, Eigen::RowMajor>(&array_F[0], fCount, 3);
        C = Eigen::Map<Eigen::MatrixXd, Eigen::RowMajor>(&array_C[0], cCount, 3);
        BE = Eigen::Map<Eigen::MatrixXi, Eigen::RowMajor>(&array_BE[0], beCount, 2);
        
        igl::boundary_conditions(V,F,C,Eigen::VectorXi(),BE,Eigen::MatrixXi(),b,bc);
        
        // compute BBW weights matrix
        igl::BBWData bbw_data;
        // only a few iterations for sake of demo
        bbw_data.active_set_params.max_iter = 8;
        bbw_data.verbosity = 2;
        
        if(!igl::bbw(V,F,b,bc,bbw_data,W))
            return EXIT_FAILURE;
        
        // Normalize weights to sum to one
        // It seems that weights does not need to sum to one in Anima2D
        //igl::normalize_col_sums(W,W);
        
        double array_W[W.rows() * W.cols()];
        Eigen::Map<Eigen::MatrixXd>(array_W, W.rows(), W.cols()) = W;
        
        for (int i = 0; i < oriVCount; i++)
            for (int j = 0; j < beCount; j++)
                ptr_W[i + oriVCount * j] = array_W[i + vCount * j];
        
        return 0;
    }
}

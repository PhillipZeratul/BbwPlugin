//
//  BbwPluginiOS.h
//  BbwPluginiOS
//
//  Created by Zeratul Phillip on 18/04/2018.
//  Copyright © 2018 Zeratul Phillip. All rights reserved.
//

// Because of Mosek complications, we don't use static library if Mosek is used.
#ifdef LIBIGL_WITH_MOSEK
#ifdef IGL_STATIC_LIBRARY
#undef IGL_STATIC_LIBRARY
#endif
#endif

#include <igl/boundary_conditions.h>
#include <igl/bbw.h>
#include <igl/normalize_row_sums.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <vector>
#include <algorithm>
#include <iostream>

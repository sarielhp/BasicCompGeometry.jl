/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
 * gdiam_test.cpp -
 *     Test code for computing diameter and tight-fitting
 * bounding box.
 *
 * Copyright 2000 Sariel Har-Peled (ssaarriieell@cs.uiuc.edu)
 *
 * * the GNU General Public License as published by the Free
 *   Software Foundation; either version 2, or (at your option)
 *   any later version.
 *
 * or
 *
 * * the GNU Lesser General Public License as published by the Free
 *   Software Foundation; either version 2.1, or (at your option)
 *   any later version.
 *
 * Code is based on the paper:
 *   A Practical Approach for Computing the Diameter of a Point-Set.
 *   Sariel Har-Peled (http://www.uiuc.edu/~sariel)
 *---------------------------------------------------------------
 * History:
 * 3/28/01 - Extended test program, so that it can load a
 *           text-file with points. Format is
 *                    [npoints]
 *                    [x_1 y_1 z_1]
 *                     ...
 *                    [x_n y_n z_n]
\*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/

#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <memory.h>
#include  <math.h>

#include  <vector>
#include  <algorithm>
#include  <random> 

#include  "mvbb.h"

/*--- Start of Code ---*/


void   test_itself( MVBB::real  * points, int  num )
{
  MVBB::GPointPair   pair;

    printf( "Computing the diameter for %d points selected "
            "uniformly from the unit cube\n", num );
    pair = MVBB::approx_diam_pair( (MVBB::real *)points, num, 0.0 );
    printf( "Diameter distance: %g\n", pair.distance );
    printf( "Points realizing the diameter\n"
            "\t(%g, %g, %g) - (%g, %g, %g)\n",
            pair.p[ 0 ], pair.p[ 1 ], pair.p[ 2 ],
            pair.q[ 0 ], pair.q[ 1 ], pair.q[ 2 ] );


    MVBB::point  * pnt_arr;
    MVBB::BBox   bb;

    pnt_arr = MVBB::convert( (MVBB::real *)points, num );

    printf( "Computing a tight-fitting bounding box of the point-set\n" );
    bb = MVBB::approx_mvbb_grid_sample( pnt_arr, num, 5, 400 );

    printf( "Resulting bounding box:\n" );
    bb.dump();

    printf( "Axis parallel bounding box\n" );
    MVBB::AABBox   bbx;
    bbx.init();
    for  ( int  ind = 0; ind < num; ind++ )
        bbx.bound( points + (ind * 3) );
    bbx.dump();
        
}


auto  rand_points( int n ) {
  auto points = std::vector<MVBB::point2d>();
  points.reserve(n); // Optimization: Reserve memory to prevent reallocations

  std::random_device rd; 
  std::mt19937 gen(rd()); 
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  for (int i = 0; i < n; ++i) {
    points.push_back( MVBB::point2d( dist(gen), dist(gen) ) );
  }
  
  return points;
}

typedef MVBB::point2d  point2d;

void   test_2d_ch()
{
  auto points = std::vector<MVBB::point2d>();
  auto ch = MVBB::vec_ptr_points2d();

  points.push_back( point2d( 0.0, 0.5 ) );
  points.push_back( point2d( 0.0, 2.5 ) );
  points.push_back( point2d( 1.0, 1.0 ) );
  points.push_back( point2d( 1.0, 1.0 ) );
  points.push_back( point2d( 0.0, 1.0 ) );
  points.push_back( point2d( 1.0, 0.0 ) );
  points.push_back( point2d( 0.0, 0.0 ) );
  points.push_back( point2d( 0.0, 0.0 ) );
  points.push_back( point2d( 0.0, 0.5 ) );
  points.push_back( point2d( 0.0, 2.5 ) );
  points.push_back( point2d( 0.0, 0.0 ) );
  points.push_back( point2d( 0.0, 0.0 ) );
  points.push_back( point2d( 0.0, 2.2498 ) );

  auto  pnt_arr = points_2d_to_pointers( points );

  convex_hull_2d( pnt_arr, ch );

  //dump( ch );
  //printf( "ch.size: %d", (int)ch.size() );
  //printf( "\n\n\n" );

  verify_convex_hull_2d( pnt_arr, ch );
  //dump( ch );
  
  //exit( -1 );

  auto pnts = rand_points( 100000 );
  auto pnts_arr = points_2d_to_pointers( pnts );
  auto ch2 = MVBB::vec_ptr_points2d();

  convex_hull_2d( pnts_arr, ch2 );
  verify_convex_hull_2d( pnts_arr, ch2 );
}



void   standard_test()
{
  test_2d_ch();

    MVBB::real  * points;
    int  num;

    num = 1000000;

    points = (MVBB::point)malloc( sizeof( MVBB::point_t ) * num );
    assert( points != NULL );

    // Pick randomly points from the unit cube */
    for  ( int  ind = 0; ind < num; ind++ ) {
        points[ ind * 3 + 0 ] = drand48();
        points[ ind * 3 + 1 ] = drand48();
        points[ ind * 3 + 2 ] = drand48();
    }

    test_itself( points, num );
}



void  read_points( FILE   * fl, MVBB::real  * points, int  points_num )
{
    int  args;
    double  x, y, z;

    for  ( int  ind = 0; ind < points_num; ind++ ) {
        args = fscanf( fl, "%lg %lg %lg\n", &x, &y, &z );
        assert( args == 3 );

        points[ ind * 3 + 0 ] = x;
        points[ ind * 3 + 1 ] = y;
        points[ ind * 3 + 2 ] = z;
    }
}


void  test_file( const char  * file_name )
{
    MVBB::real  * points;
    FILE   * fl;
    int  args, points_num;

    fl = fopen( file_name, "rt" );
    if  ( fl == NULL ) {
        printf( "Unable to open file: [%s]\n", file_name );
        exit( -1 );
    }
    args = fscanf( fl, "%d\n", &points_num );
    assert( ( args > 0 )  &&  ( points_num > 0 ) );

    points = (MVBB::point)malloc( sizeof( MVBB::point_t ) * points_num );
    assert( points != NULL );

    read_points( fl, points, points_num );
    fclose( fl );

    test_itself( points, points_num );
}


int  main( int  argc, char  ** argv )
{
    if  ( argc == 1 ) {
        standard_test();
        return  0;
    }

    for  ( int  ind = 1; ind < argc; ind++ )
        test_file( argv[ ind ] );

    return  0;
}

/* mvbb_test.C - End of File ------------------------------------------*/

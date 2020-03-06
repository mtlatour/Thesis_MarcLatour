/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "codedFvOptionTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "fvMatrix.H"

//{{{ begin codeInclude
#line 30 "/home/palak/OpenFOAM/palak-6/run/Z4curvedlookuptable/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
#include </usr/include/CGAL/Exact_predicates_inexact_constructions_kernel.h>
				#include </usr/include/CGAL/Delaunay_triangulation_2.h>
				#include </usr/include/CGAL/Interpolation_traits_2.h>
				#include </usr/include/CGAL/natural_neighbor_coordinates_2.h>
				#include </usr/include/CGAL/interpolation_functions.h>
				#include "IFstream.H"
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 90e2d38f924a8f1f59cf06935fc5b1741620a6f6
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void bodyForce_90e2d38f924a8f1f59cf06935fc5b1741620a6f6(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//makeRemovablePatchTypeField
//(
//    fvPatch,
//    bodyForceFvOptionvectorSource
//);
defineTypeNameAndDebug(bodyForceFvOptionvectorSource, 0);
addRemovableToRunTimeSelectionTable
(
    option,
    bodyForceFvOptionvectorSource,
    dictionary
);


const char* const bodyForceFvOptionvectorSource::SHA1sum =
    "90e2d38f924a8f1f59cf06935fc5b1741620a6f6";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

bodyForceFvOptionvectorSource::
bodyForceFvOptionvectorSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh)
{
    if (false)
    {
        Info<<"construct bodyForce sha1: 90e2d38f924a8f1f59cf06935fc5b1741620a6f6"
            " from components\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

bodyForceFvOptionvectorSource::
~bodyForceFvOptionvectorSource()
{
    if (false)
    {
        Info<<"destroy bodyForce sha1: 90e2d38f924a8f1f59cf06935fc5b1741620a6f6\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void bodyForceFvOptionvectorSource::correct
(
    GeometricField<vector, fvPatchField, volMesh>& fld
)
{
    if (false)
    {
        Info<<"bodyForceFvOptionvectorSource::correct()\n";
    }

//{{{ begin code
    #line 39 "/home/palak/OpenFOAM/palak-6/run/Z4curvedlookuptable/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
Pout<< "**codeCorrect**" << endl;
//}}} end code
}


void bodyForceFvOptionvectorSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (false)
    {
        Info<<"bodyForceFvOptionvectorSource::addSup()\n";
    }

//{{{ begin code
    #line 43 "/home/palak/OpenFOAM/palak-6/run/Z4curvedlookuptable/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
//Known constant values: B, R, omega, theta
				const scalar B = 1/2.5;
				const scalar R = 1;
				const scalar omega = 0;
				scalar PI = constant::mathematical::pi;
				scalar OMEGAR = omega*R;
				//Initializing all fields
				
				int s = cells_.size();

				scalarField slope(cells_.size());
				scalarField angle(cells_.size());
				scalarField x(s), y(s), z(s), RADIUS(s), THETA(s), NX(s), NY(s), NZ(s), NR(s), NTH(s), WX(s), WY(s), WZ(s), 
							WMAG(s), WDOTN(s), WNX(s), WNY(s), WNZ(s), DEVLOC(s), WTX(s), WTY(s), WTZ(s), WTMAG(s), TX(s), TY(s), TZ(s), FNX(s), FNY(s), FNZ(s),
							FTX(s), FTY(s), FTZ(s), MOMSRCX(s), MOMSRCY(s), MOMSRCZ(s);
				
				volVectorField bodyForce
				(
					IOobject
					(
						name_ + ":bodyForce",
						mesh_.time().timeName(),
						mesh_
					),
					mesh_,
					dimensionedVector
					(
						"zero",
						eqn.dimensions()/dimVolume,
						Zero
					)
				);
				
				
				const vectorField& U = eqn.psi();
				const vectorField& CC = mesh_.C(); //cell center 
					
				//SETTING UP LINEAR INTERPOLATION USING CGAL FUNCTIONS
		
				typedef CGAL::Exact_predicates_inexact_constructions_kernel 	K;
				typedef CGAL::Delaunay_triangulation_2<K>                   	Delaunay_triangulation;
				typedef K::FT                                               	Coord_type;
				typedef K::Point_2                                          	Point;
				typedef std::map<Point, Coord_type, K::Less_xy_2>         		Coord_map;
				typedef CGAL::Data_access<Coord_map>                      		Value_access;
				
				Delaunay_triangulation Tnx, Tny;										//Holds the points
				
				Coord_map value_nx, value_ny;                                       //Holds the known points and their known corresponding values
				
					List<vector> nx_data, ny_data;
					IFstream nx("nx_data");
					nx  >> nx_data;
					IFstream ny("ny_data");
					ny  >> ny_data;

				forAll(nx_data, a)
				{
					K::Point_2 px(nx_data[a][0],nx_data[a][1]);
					Tnx.insert(px);
					value_nx.insert(std::make_pair(px,nx_data[a][2]));
				}

				forAll(ny_data, b)
				{
					K::Point_2 py(ny_data[b][0],ny_data[b][1]);
					Tny.insert(py);
					value_ny.insert(std::make_pair(py,ny_data[b][2]));
				}
					
				
				forAll(cells_, m)
				{
				
					x[m] = CC[m].x();
					y[m] = CC[m].y();
					K::Point_2 p(x[m], y[m]);
					std::vector<std::pair<Point, Coord_type> > coords_nx, coords_ny;
					Coord_type norm_nx = CGAL::natural_neighbor_coordinates_2(Tnx, p, std::back_inserter(coords_nx)).second;
					Coord_type norm_ny = CGAL::natural_neighbor_coordinates_2(Tny, p, std::back_inserter(coords_ny)).second;
					Coord_type res_nx =  CGAL::linear_interpolation(coords_nx.begin(), coords_nx.end(), norm_nx, Value_access(value_nx));
					Coord_type res_ny =  CGAL::linear_interpolation(coords_ny.begin(), coords_ny.end(), norm_ny, Value_access(value_ny));
					NX[m]=res_nx;
					NY[m]=res_ny;
				}
				
				forAll(cells_, i)
				{
					if((CC[i].x()>=2.0) && (CC[i].x()<=4.5))

					{ 
					
					WX[i]= U[i].x();
					WY[i]= U[i].y() + OMEGAR;
					WZ[i]= 0;
					WMAG[i] = sqrt(WX[i]*WX[i] + WY[i]*WY[i] + WZ[i]*WZ[i]); 
					WDOTN[i] = WX[i]*NX[i] + WY[i]*NY[i] + WZ[i]*NZ[i];
					DEVLOC[i] = asin(WDOTN[i]/WMAG[i]);
					FNX[i] = -PI*B*DEVLOC[i]*WMAG[i]*WMAG[i]*sin(DEVLOC[i]);
					FNY[i] = -PI*B*DEVLOC[i]*WMAG[i]*WMAG[i]*cos(DEVLOC[i]);
					FNZ[i] = 0;
					MOMSRCX[i] = (FNX[i]);
					MOMSRCY[i] = (FNY[i]);
					MOMSRCZ[i] = (FNZ[i]);
					bodyForce[i] = vector(MOMSRCX[i], MOMSRCY[i], MOMSRCZ[i]);
					
					}
					
					else
					
					{
					
						MOMSRCX[i] = 0;
						MOMSRCY[i] = 0;
						MOMSRCZ[i] = 0;
					
						bodyForce[i] = vector(0, 0, 0);
					
					}
					
				// 	if(mesh_.time().outputTime())
				// 	{
				// 		std::cout << " At ("	<< CC[i].x() << ","<< CC[i].y() << ","<< CC[i].z() << ") " << 
				// 		" Deviation = " << DEVLOC[i] <<
				// 		" Unit vector = (" << NX[i] <<  ","	<< NY[i] << ")" <<
				// 		" F = (" << MOMSRCX[i] << ","	<< MOMSRCY[i] << ","<< MOMSRCZ[i] << ") " << std::endl;
				// 		bodyForce.write();
				// 	}
				 }
			
					
          	        eqn += bodyForce;
//}}} end code
}


void bodyForceFvOptionvectorSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (false)
    {
        Info<<"bodyForceFvOptionvectorSource::addSup()\n";
    }

//{{{ begin code
    #line 43 "/home/palak/OpenFOAM/palak-6/run/Z4curvedlookuptable/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
//Known constant values: B, R, omega, theta
				const scalar B = 1/2.5;
				const scalar R = 1;
				const scalar omega = 0;
				scalar PI = constant::mathematical::pi;
				scalar OMEGAR = omega*R;
				//Initializing all fields
				
				int s = cells_.size();

				scalarField slope(cells_.size());
				scalarField angle(cells_.size());
				scalarField x(s), y(s), z(s), RADIUS(s), THETA(s), NX(s), NY(s), NZ(s), NR(s), NTH(s), WX(s), WY(s), WZ(s), 
							WMAG(s), WDOTN(s), WNX(s), WNY(s), WNZ(s), DEVLOC(s), WTX(s), WTY(s), WTZ(s), WTMAG(s), TX(s), TY(s), TZ(s), FNX(s), FNY(s), FNZ(s),
							FTX(s), FTY(s), FTZ(s), MOMSRCX(s), MOMSRCY(s), MOMSRCZ(s);
				
				volVectorField bodyForce
				(
					IOobject
					(
						name_ + ":bodyForce",
						mesh_.time().timeName(),
						mesh_
					),
					mesh_,
					dimensionedVector
					(
						"zero",
						eqn.dimensions()/dimVolume,
						Zero
					)
				);
				
				
				const vectorField& U = eqn.psi();
				const vectorField& CC = mesh_.C(); //cell center 
					
				//SETTING UP LINEAR INTERPOLATION USING CGAL FUNCTIONS
		
				typedef CGAL::Exact_predicates_inexact_constructions_kernel 	K;
				typedef CGAL::Delaunay_triangulation_2<K>                   	Delaunay_triangulation;
				typedef K::FT                                               	Coord_type;
				typedef K::Point_2                                          	Point;
				typedef std::map<Point, Coord_type, K::Less_xy_2>         		Coord_map;
				typedef CGAL::Data_access<Coord_map>                      		Value_access;
				
				Delaunay_triangulation Tnx, Tny;										//Holds the points
				
				Coord_map value_nx, value_ny;                                       //Holds the known points and their known corresponding values
				
					List<vector> nx_data, ny_data;
					IFstream nx("nx_data");
					nx  >> nx_data;
					IFstream ny("ny_data");
					ny  >> ny_data;

				forAll(nx_data, a)
				{
					K::Point_2 px(nx_data[a][0],nx_data[a][1]);
					Tnx.insert(px);
					value_nx.insert(std::make_pair(px,nx_data[a][2]));
				}

				forAll(ny_data, b)
				{
					K::Point_2 py(ny_data[b][0],ny_data[b][1]);
					Tny.insert(py);
					value_ny.insert(std::make_pair(py,ny_data[b][2]));
				}
					
				
				forAll(cells_, m)
				{
				
					x[m] = CC[m].x();
					y[m] = CC[m].y();
					K::Point_2 p(x[m], y[m]);
					std::vector<std::pair<Point, Coord_type> > coords_nx, coords_ny;
					Coord_type norm_nx = CGAL::natural_neighbor_coordinates_2(Tnx, p, std::back_inserter(coords_nx)).second;
					Coord_type norm_ny = CGAL::natural_neighbor_coordinates_2(Tny, p, std::back_inserter(coords_ny)).second;
					Coord_type res_nx =  CGAL::linear_interpolation(coords_nx.begin(), coords_nx.end(), norm_nx, Value_access(value_nx));
					Coord_type res_ny =  CGAL::linear_interpolation(coords_ny.begin(), coords_ny.end(), norm_ny, Value_access(value_ny));
					NX[m]=res_nx;
					NY[m]=res_ny;
				}
				
				forAll(cells_, i)
				{
					if((CC[i].x()>=2.0) && (CC[i].x()<=4.5))

					{ 
					
					WX[i]= U[i].x();
					WY[i]= U[i].y() + OMEGAR;
					WZ[i]= 0;
					WMAG[i] = sqrt(WX[i]*WX[i] + WY[i]*WY[i] + WZ[i]*WZ[i]); 
					WDOTN[i] = WX[i]*NX[i] + WY[i]*NY[i] + WZ[i]*NZ[i];
					DEVLOC[i] = asin(WDOTN[i]/WMAG[i]);
					FNX[i] = -PI*B*DEVLOC[i]*WMAG[i]*WMAG[i]*sin(DEVLOC[i]);
					FNY[i] = -PI*B*DEVLOC[i]*WMAG[i]*WMAG[i]*cos(DEVLOC[i]);
					FNZ[i] = 0;
					MOMSRCX[i] = (FNX[i]);
					MOMSRCY[i] = (FNY[i]);
					MOMSRCZ[i] = (FNZ[i]);
					bodyForce[i] = vector(MOMSRCX[i], MOMSRCY[i], MOMSRCZ[i]);
					
					}
					
					else
					
					{
					
						MOMSRCX[i] = 0;
						MOMSRCY[i] = 0;
						MOMSRCZ[i] = 0;
					
						bodyForce[i] = vector(0, 0, 0);
					
					}
					
				// 	if(mesh_.time().outputTime())
				// 	{
				// 		std::cout << " At ("	<< CC[i].x() << ","<< CC[i].y() << ","<< CC[i].z() << ") " << 
				// 		" Deviation = " << DEVLOC[i] <<
				// 		" Unit vector = (" << NX[i] <<  ","	<< NY[i] << ")" <<
				// 		" F = (" << MOMSRCX[i] << ","	<< MOMSRCY[i] << ","<< MOMSRCZ[i] << ") " << std::endl;
				// 		bodyForce.write();
				// 	}
				 }
			
					
          	        eqn += bodyForce;
//}}} end code
}


void bodyForceFvOptionvectorSource::constrain
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (false)
    {
        Info<<"bodyForceFvOptionvectorSource::constrain()\n";
    }

//{{{ begin code
    #line 179 "/home/palak/OpenFOAM/palak-6/run/Z4curvedlookuptable/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
Pout<< "**codeSetValue**" << endl;
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

} // End namespace fv
// ************************************************************************* //


/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

Description
    Template for use with dynamic code generation of a source.
    The hook functions take the following arguments:

    codeCorrect
    (
        GeometricField<Type, fvPatchField, volMesh>& fld
    )

    codeAddSup
    (
        fvMatrix<Type}>& eqn,
        const label fieldi
    )

    constrain
    (
        fvMatrix<Type}>& eqn,
        const label fieldi
    )

    where :
        fld is the field in fieldNames
        eqn is the fvMatrix

    energySource
    {
        type            scalarCodedSource;

        active          true;
        selectionMode   all;

        scalarCodedSourceCoeffs
        {
            fieldNames      (h);
            name    sourceTime;

            codeInclude
            #{

            #};

            codeCorrect
            #{
                Pout<< "**codeCorrect**" << endl;
            #};

            codeAddSup
            #{
                const Time& time = mesh().time();
                const scalarField& V = mesh_.V();
                scalarField& heSource = eqn.source();
                heSource -= 0.1*sqr(time.value())*V;
            #};

            codeSetValue
            #{
                Pout<< "**codeSetValue**" << endl;
            #};

            // Dummy entry. Make dependent on above to trigger recompilation
            code
            #{
                #line 30 "/home/palak/OpenFOAM/palak-5.0/run/Z4curvedlookuptable/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
#include </usr/include/CGAL/Exact_predicates_inexact_constructions_kernel.h>
				#include </usr/include/CGAL/Delaunay_triangulation_2.h>
				#include </usr/include/CGAL/Interpolation_traits_2.h>
				#include </usr/include/CGAL/natural_neighbor_coordinates_2.h>
				#include </usr/include/CGAL/interpolation_functions.h>
				#include "IFstream.H"
                #line 39 "/home/palak/OpenFOAM/palak-5.0/run/Z4curvedlookuptable/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
Pout<< "**codeCorrect**" << endl;
                #line 43 "/home/palak/OpenFOAM/palak-5.0/run/Z4curvedlookuptable/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
//Known constant values: B, R, omega, theta
				const scalar B = 5;
				const scalar R = 1;
				const scalar omega = 0;
				const scalar theta = 10;
				scalar PI = constant::mathematical::pi;
				scalar OMEGAR = omega*R;
				
				//Initializing all fields
				
				scalarField x_coord(cells_.size());
				scalarField y_coord(cells_.size());
				scalarField NX(cells_.size());
				scalarField NY(cells_.size());
				scalarField NZ(cells_.size());
				scalarField slope(cells_.size());
				scalarField angle(cells_.size());
				scalarField WX(cells_.size());
				scalarField WY(cells_.size());
				scalarField WZ(cells_.size());
				scalarField WMAG(cells_.size());
				scalarField WDOTN(cells_.size());
				scalarField WNX(cells_.size());
				scalarField WNY(cells_.size());
				scalarField WNZ(cells_.size());
				scalarField DEVLOC(cells_.size());
				scalarField FNX(cells_.size());
				scalarField FNY(cells_.size());
				scalarField FNZ(cells_.size());
				scalarField MOMSRCX(cells_.size());
				scalarField MOMSRCY(cells_.size());
				scalarField MOMSRCZ(cells_.size());
				
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
				
				Delaunay_triangulation Tnx;										//Holds the points
				Delaunay_triangulation Tny;
				
				Coord_map value_function;                                       //Holds the known points and their known corresponding values
				Coord_map value_ny;  
				
					List<vector> nx_data;
					IFstream nx("nx_data");
					nx  >> nx_data;
					
					List<vector> ny_data;
					IFstream ny("ny_data");
					ny  >> ny_data;

					if (nx_data.size() > 0)
					{
						scalarField Xnx(nx_data.size());
						scalarField Ynx(nx_data.size());
						scalarField nx(nx_data.size());

						forAll(nx_data, j)
						{
							Xnx[j] = nx_data[j][0];
							Ynx[j] = nx_data[j][1];
							 nx[j] = nx_data[j][2];
							
							K::Point_2 p(Xnx[j],Ynx[j]);
							Tnx.insert(p);
							value_function.insert(std::make_pair(p, nx[j]));
							//Outputting data on each timestep		
							//std::cout << "x " << Xnx[j] << " y "	<< Ynx[j] << " nx "<< nx[j] << std::endl; 
						
						}
					}
					
				
				forAll(cells_, m)
				{
				
					x_coord[m] = CC[m].x();
					y_coord[m] = CC[m].y();
					
					K::Point_2 p(x_coord[m], y_coord[m]);
					std::vector<std::pair<Point, Coord_type> > coords;
					Coord_type norm = CGAL::natural_neighbor_coordinates_2(Tnx, p, std::back_inserter(coords)).second;
					Coord_type res =  CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(value_function));
					NX[m]=res;
				}
					
					
					
					
					if (ny_data.size() > 0)
					{
						scalarField Xny(ny_data.size());
						scalarField Yny(ny_data.size());
						scalarField  ny(ny_data.size());

						forAll(ny_data, k)
						{
							Xny[k] = ny_data[k][0];
							Yny[k] = ny_data[k][1];
							 ny[k] = ny_data[k][2];
							
							K::Point_2 p(Xny[k],Yny[k]);
							Tny.insert(p);
							value_ny.insert(std::make_pair(p, ny[k]));
						
						}
					}
					
				
				forAll(cells_, l)
				{
					
					////coordinate computation
					K::Point_2 p(x_coord[l], y_coord[l]);
					std::vector<std::pair<Point, Coord_type> > coords2;
					Coord_type normy = CGAL::natural_neighbor_coordinates_2(Tny, p, std::back_inserter(coords2)).second;
					Coord_type resy =  CGAL::linear_interpolation(coords2.begin(), coords2.end(), normy, Value_access(value_ny));
					NY[l]=resy;
					//std::cout << "Interpolated unit vector at (" << x_coord[l] << "," << y_coord[l] << ") is" << resy << std::endl;
				}
				
				
				
				forAll(cells_, i)
				{
					if((CC[i].x()>=0.5) && (CC[i].x()<=1.5))

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
					
					if(mesh_.time().outputTime())
					{
						std::cout << " At ("	<< CC[i].x() << ","<< CC[i].y() << ","<< CC[i].z() << ") " << 
						" Deviation = " << DEVLOC[i] <<
						" Unit vector = (" << NX[i] <<  ","	<< NY[i] << ")" <<
						" F = (" << MOMSRCX[i] << ","	<< MOMSRCY[i] << ","<< MOMSRCZ[i] << ") " << std::endl;
						bodyForce.write();
					}
				}
			
					
          	        eqn += bodyForce;
                #line 242 "/home/palak/OpenFOAM/palak-5.0/run/Z4curvedlookuptable/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
Pout<< "**codeSetValue**" << endl;
            #};
        }

        sourceTimeCoeffs
        {
            // Dummy entry
        }
    }

SourceFiles
    codedFvOptionTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef codedFvOptionTemplate_H
#define codedFvOptionTemplate_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                         A templated CodedFvOption
\*---------------------------------------------------------------------------*/

class bodyForceFvOptionvectorSource
:
    public cellSetOption
{
public:

    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("bodyForce");


    // Constructors

        //- Construct from patch and internal field
        bodyForceFvOptionvectorSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~bodyForceFvOptionvectorSource();


    // Member functions

            //- Correct field
            virtual void correct
            (
                GeometricField<vector, fvPatchField, volMesh>&
            );

            //- Explicit and implicit matrix contributions
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Explicit and implicit matrix contributions for compressible
            //  equations
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Set value
            virtual void constrain
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


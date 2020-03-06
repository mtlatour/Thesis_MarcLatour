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
#line 30 "/home/palak/OpenFOAM/palak-6/run/Z3curvedanalytical/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
#include <IFstream.H>
            	#include <OFstream.H>
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
    // SHA1 = b70d5b00759fcdda7cfbae037b9d54fe9bba9fd8
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void bodyForce_b70d5b00759fcdda7cfbae037b9d54fe9bba9fd8(bool load)
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
    "b70d5b00759fcdda7cfbae037b9d54fe9bba9fd8";


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
        Info<<"construct bodyForce sha1: b70d5b00759fcdda7cfbae037b9d54fe9bba9fd8"
            " from components\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

bodyForceFvOptionvectorSource::
~bodyForceFvOptionvectorSource()
{
    if (false)
    {
        Info<<"destroy bodyForce sha1: b70d5b00759fcdda7cfbae037b9d54fe9bba9fd8\n";
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
    #line 36 "/home/palak/OpenFOAM/palak-6/run/Z3curvedanalytical/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
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
    #line 40 "/home/palak/OpenFOAM/palak-6/run/Z3curvedanalytical/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
//Known constant values: B, R, omega, theta
				const scalar B = 0.4;
				const scalar R = 4.75;
				const scalar omega = (1000/4.75);
				scalar PI = constant::mathematical::pi;
				scalar OMEGAR = omega*R;
				
				//Initializing all fields
				
int s = cells_.size();
				scalarField x(s), y(s), z(s), RADIUS(s), THETA(s), NX(s), NY(s), NZ(s), NR(s), NTH(s), WX(s), WY(s), WZ(s), 
							WMAG(s), WDOTN(s), WNX(s), WNY(s), WNZ(s), DEVLOC(s), WTX(s), WTY(s), WTZ(s), WTMAG(s), TX(s), TY(s), TZ(s), 
							FNX(s), FNY(s), FNZ(s), FTX(s), FTY(s), FTZ(s), MOMSRCX(s), MOMSRCY(s), MOMSRCZ(s), N(s), angle(s), x_coord(s), y_coord(s);
				
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
				
								volScalarField slope
				(
					IOobject
					(
						name_ + ":slope",
						mesh_.time().timeName(),
						mesh_
					),
					mesh_,
					dimensionedScalar
					(
						"zero",
						eqn.dimensions()/dimVolume,
						Zero
					)
				);
				
				volScalarField fangle
				(
					IOobject
					(
						name_ + ":fangle",
						mesh_.time().timeName(),
						mesh_
					),
					mesh_,
					dimensionedScalar
					(
						"zero",
						eqn.dimensions()/dimVolume,
						Zero
					)
				);
				
				volScalarField deviation
				(
					IOobject
					(
						name_ + ":deviation",
						mesh_.time().timeName(),
						mesh_
					),
					mesh_,
					dimensionedScalar
					(
						"zero",
						eqn.dimensions()/dimVolume,
						Zero
					)
				);
				
				const vectorField& U = eqn.psi();
				const vectorField& CC = mesh_.C(); //cell center 
				

					


				OFstream file("file");
  				file << "x r nx nth nr \n";


				forAll(cells_, j)
				{
				
						x_coord[j] = CC[j].x();
						y_coord[j] = CC[j].y();
						
						if((x_coord[j]>=2) && (x_coord[j]<=4.5))
						{

						NX[j] =   0.069328*x_coord[j]+0.330814;
						NY[j] =  -sqrt(1 - NX[j]*NX[j]);
						NZ[j] =  0;
						N[j] =  1/NY[j];
						
						file << x_coord[j] << "\t" << y_coord[j] << "\t" << NX[j] << "\t" << NY[j] << "\t" << endl;
						
						}
				}


				forAll(cells_, i)
				{
					if((CC[i].x()>=2) && (CC[i].x()<=4.5))

					{ 
					
					WX[i]= U[i].x();
					WY[i]= U[i].y() + OMEGAR;
					WZ[i]= 0;
					WMAG[i] = sqrt(WX[i]*WX[i] + WY[i]*WY[i] + WZ[i]*WZ[i]); 
					WDOTN[i] = WX[i]*NX[i] + WY[i]*NY[i] + WZ[i]*NZ[i];

					WNX[i] = WDOTN[i]*NX[i];
						WNY[i] = WDOTN[i]*NY[i];
						WNZ[i] = WDOTN[i]*NZ[i];


					DEVLOC[i] = asin(WDOTN[i]/WMAG[i]);
											
						WTX[i] = WX[i] - WNX[i];
						WTY[i] = WY[i] - WNY[i];
						WTZ[i] = WZ[i] - WNZ[i];
						WTMAG[i] = sqrt(WTX[i]*WTX[i] + WTY[i]*WTY[i] + WTZ[i]*WTZ[i]);
						
						TX[i] = WTX[i]/WTMAG[i];
						TY[i] = WTY[i]/WTMAG[i];
						TZ[i] = WTZ[i]/WTMAG[i];
						
						FNX[i] = -2*PI*DEVLOC[i]*cos(DEVLOC[i])*WMAG[i]*WMAG[i]*B*NX[i];
						FNY[i] = -2*PI*DEVLOC[i]*cos(DEVLOC[i])*WMAG[i]*WMAG[i]*B*NY[i];
						FNZ[i] = -2*PI*DEVLOC[i]*cos(DEVLOC[i])*WMAG[i]*WMAG[i]*B*NZ[i];

						FTX[i] = 2*PI*DEVLOC[i]*sin(DEVLOC[i])*WMAG[i]*WMAG[i]*B*TX[i];
						FTY[i] = 2*PI*DEVLOC[i]*sin(DEVLOC[i])*WMAG[i]*WMAG[i]*B*TY[i];
						FTZ[i] = 2*PI*DEVLOC[i]*sin(DEVLOC[i])*WMAG[i]*WMAG[i]*B*TZ[i];
						
						MOMSRCX[i] = (FTX[i] + FNX[i]);
						MOMSRCY[i] = (FTY[i] + FNY[i]);
						MOMSRCZ[i] = (FTZ[i] + FNZ[i]);

					// 	FNX[i] = -PI*B*DEVLOC[i]*WMAG[i]*WMAG[i]*cos(DEVLOC[i]);
					// FNY[i] = -PI*B*DEVLOC[i]*WMAG[i]*WMAG[i]*sin(DEVLOC[i]);
					// FNZ[i] = 0;
					// MOMSRCX[i] = (FNX[i]);
					// MOMSRCY[i] = (FNY[i]);
					// MOMSRCZ[i] = (FNZ[i]);
					bodyForce[i] = vector(MOMSRCX[i], MOMSRCY[i], MOMSRCZ[i]);
					
					
					}
					
					else
					
					{
					
						MOMSRCX[i] = 0;
						MOMSRCY[i] = 0;
						MOMSRCZ[i] = 0;
					
						bodyForce[i] = vector(0, 0, 0);
					
					}
					
					//if(mesh_.time().outputTime())
					//{
					//	std::cout << " At ("	<< CC[i].x() << ","<< CC[i].y() << ","<< CC[i].z() << ") " << 
					//	" Deviation = " << DEVLOC[i] <<
					//	" Unit vector = (" << NX[i] <<  ","	<< NY[i] << ")" <<
					//	" F = (" << MOMSRCX[i] << ","	<< MOMSRCY[i] << ","<< MOMSRCZ[i] << ") " << std::endl;
					//	bodyForce.write();
					//}
				}
			
					
          	        eqn += bodyForce;
					
					if(mesh_.time().outputTime())
					{
						bodyForce.write();
						slope.write();
						fangle.write();
						deviation.write();
					}
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
    #line 40 "/home/palak/OpenFOAM/palak-6/run/Z3curvedanalytical/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
//Known constant values: B, R, omega, theta
				const scalar B = 0.4;
				const scalar R = 4.75;
				const scalar omega = (1000/4.75);
				scalar PI = constant::mathematical::pi;
				scalar OMEGAR = omega*R;
				
				//Initializing all fields
				
int s = cells_.size();
				scalarField x(s), y(s), z(s), RADIUS(s), THETA(s), NX(s), NY(s), NZ(s), NR(s), NTH(s), WX(s), WY(s), WZ(s), 
							WMAG(s), WDOTN(s), WNX(s), WNY(s), WNZ(s), DEVLOC(s), WTX(s), WTY(s), WTZ(s), WTMAG(s), TX(s), TY(s), TZ(s), 
							FNX(s), FNY(s), FNZ(s), FTX(s), FTY(s), FTZ(s), MOMSRCX(s), MOMSRCY(s), MOMSRCZ(s), N(s), angle(s), x_coord(s), y_coord(s);
				
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
				
								volScalarField slope
				(
					IOobject
					(
						name_ + ":slope",
						mesh_.time().timeName(),
						mesh_
					),
					mesh_,
					dimensionedScalar
					(
						"zero",
						eqn.dimensions()/dimVolume,
						Zero
					)
				);
				
				volScalarField fangle
				(
					IOobject
					(
						name_ + ":fangle",
						mesh_.time().timeName(),
						mesh_
					),
					mesh_,
					dimensionedScalar
					(
						"zero",
						eqn.dimensions()/dimVolume,
						Zero
					)
				);
				
				volScalarField deviation
				(
					IOobject
					(
						name_ + ":deviation",
						mesh_.time().timeName(),
						mesh_
					),
					mesh_,
					dimensionedScalar
					(
						"zero",
						eqn.dimensions()/dimVolume,
						Zero
					)
				);
				
				const vectorField& U = eqn.psi();
				const vectorField& CC = mesh_.C(); //cell center 
				

					


				OFstream file("file");
  				file << "x r nx nth nr \n";


				forAll(cells_, j)
				{
				
						x_coord[j] = CC[j].x();
						y_coord[j] = CC[j].y();
						
						if((x_coord[j]>=2) && (x_coord[j]<=4.5))
						{

						NX[j] =   0.069328*x_coord[j]+0.330814;
						NY[j] =  -sqrt(1 - NX[j]*NX[j]);
						NZ[j] =  0;
						N[j] =  1/NY[j];
						
						file << x_coord[j] << "\t" << y_coord[j] << "\t" << NX[j] << "\t" << NY[j] << "\t" << endl;
						
						}
				}


				forAll(cells_, i)
				{
					if((CC[i].x()>=2) && (CC[i].x()<=4.5))

					{ 
					
					WX[i]= U[i].x();
					WY[i]= U[i].y() + OMEGAR;
					WZ[i]= 0;
					WMAG[i] = sqrt(WX[i]*WX[i] + WY[i]*WY[i] + WZ[i]*WZ[i]); 
					WDOTN[i] = WX[i]*NX[i] + WY[i]*NY[i] + WZ[i]*NZ[i];

					WNX[i] = WDOTN[i]*NX[i];
						WNY[i] = WDOTN[i]*NY[i];
						WNZ[i] = WDOTN[i]*NZ[i];


					DEVLOC[i] = asin(WDOTN[i]/WMAG[i]);
											
						WTX[i] = WX[i] - WNX[i];
						WTY[i] = WY[i] - WNY[i];
						WTZ[i] = WZ[i] - WNZ[i];
						WTMAG[i] = sqrt(WTX[i]*WTX[i] + WTY[i]*WTY[i] + WTZ[i]*WTZ[i]);
						
						TX[i] = WTX[i]/WTMAG[i];
						TY[i] = WTY[i]/WTMAG[i];
						TZ[i] = WTZ[i]/WTMAG[i];
						
						FNX[i] = -2*PI*DEVLOC[i]*cos(DEVLOC[i])*WMAG[i]*WMAG[i]*B*NX[i];
						FNY[i] = -2*PI*DEVLOC[i]*cos(DEVLOC[i])*WMAG[i]*WMAG[i]*B*NY[i];
						FNZ[i] = -2*PI*DEVLOC[i]*cos(DEVLOC[i])*WMAG[i]*WMAG[i]*B*NZ[i];

						FTX[i] = 2*PI*DEVLOC[i]*sin(DEVLOC[i])*WMAG[i]*WMAG[i]*B*TX[i];
						FTY[i] = 2*PI*DEVLOC[i]*sin(DEVLOC[i])*WMAG[i]*WMAG[i]*B*TY[i];
						FTZ[i] = 2*PI*DEVLOC[i]*sin(DEVLOC[i])*WMAG[i]*WMAG[i]*B*TZ[i];
						
						MOMSRCX[i] = (FTX[i] + FNX[i]);
						MOMSRCY[i] = (FTY[i] + FNY[i]);
						MOMSRCZ[i] = (FTZ[i] + FNZ[i]);

					// 	FNX[i] = -PI*B*DEVLOC[i]*WMAG[i]*WMAG[i]*cos(DEVLOC[i]);
					// FNY[i] = -PI*B*DEVLOC[i]*WMAG[i]*WMAG[i]*sin(DEVLOC[i]);
					// FNZ[i] = 0;
					// MOMSRCX[i] = (FNX[i]);
					// MOMSRCY[i] = (FNY[i]);
					// MOMSRCZ[i] = (FNZ[i]);
					bodyForce[i] = vector(MOMSRCX[i], MOMSRCY[i], MOMSRCZ[i]);
					
					
					}
					
					else
					
					{
					
						MOMSRCX[i] = 0;
						MOMSRCY[i] = 0;
						MOMSRCZ[i] = 0;
					
						bodyForce[i] = vector(0, 0, 0);
					
					}
					
					//if(mesh_.time().outputTime())
					//{
					//	std::cout << " At ("	<< CC[i].x() << ","<< CC[i].y() << ","<< CC[i].z() << ") " << 
					//	" Deviation = " << DEVLOC[i] <<
					//	" Unit vector = (" << NX[i] <<  ","	<< NY[i] << ")" <<
					//	" F = (" << MOMSRCX[i] << ","	<< MOMSRCY[i] << ","<< MOMSRCZ[i] << ") " << std::endl;
					//	bodyForce.write();
					//}
				}
			
					
          	        eqn += bodyForce;
					
					if(mesh_.time().outputTime())
					{
						bodyForce.write();
						slope.write();
						fangle.write();
						deviation.write();
					}
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
    #line 241 "/home/palak/OpenFOAM/palak-6/run/Z3curvedanalytical/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
Pout<< "**codeSetValue**" << endl;
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

} // End namespace fv
// ************************************************************************* //


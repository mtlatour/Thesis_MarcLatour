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
#line 30 "/home/palak/OpenFOAM/palak-6/run/Z1flatanalytical/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"

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
    // SHA1 = c4d8b43b18a37e1160aa3b34285a1fee1b55d330
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void bodyForce_c4d8b43b18a37e1160aa3b34285a1fee1b55d330(bool load)
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
    "c4d8b43b18a37e1160aa3b34285a1fee1b55d330";


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
        Info<<"construct bodyForce sha1: c4d8b43b18a37e1160aa3b34285a1fee1b55d330"
            " from components\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

bodyForceFvOptionvectorSource::
~bodyForceFvOptionvectorSource()
{
    if (false)
    {
        Info<<"destroy bodyForce sha1: c4d8b43b18a37e1160aa3b34285a1fee1b55d330\n";
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
    #line 33 "/home/palak/OpenFOAM/palak-6/run/Z1flatanalytical/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
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
    #line 37 "/home/palak/OpenFOAM/palak-6/run/Z1flatanalytical/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
//Known constant values: B, R, omega, theta
				const scalar B = 1;
				const scalar R = 1;
				const scalar omega = 1000;
				const scalar theta = 40;
				
				//Converting theta to radians
				scalar PI = constant::mathematical::pi;
				scalar radTHETA = (theta*PI/180.0);
				
				//Calculating OMEGAR, NX, NY, NZ
				scalar OMEGAR = omega*R;
				const scalar NX = sin(radTHETA);
				const scalar NY = cos(radTHETA);
				const scalar NZ = 0;
				
				//Initializing all fields
				
				scalarField x_coord(cells_.size());
				scalarField y_coord(cells_.size());
				
				scalarField WX(cells_.size());
				scalarField WY(cells_.size());
				scalarField WZ(cells_.size());
				scalarField WMAG(cells_.size());
				
				scalarField WDOTN(cells_.size());
				
				scalarField WNX(cells_.size());
				scalarField WNY(cells_.size());
				scalarField WNZ(cells_.size());
				scalarField DEVLOC(cells_.size());
				scalarField angle(cells_.size());
				
				scalarField WTX(cells_.size());
				scalarField WTY(cells_.size());
				scalarField WTZ(cells_.size());
				scalarField WTMAG(cells_.size());
				
				scalarField TX(cells_.size());
				scalarField TY(cells_.size());
				scalarField TZ(cells_.size());
				
				scalarField FNX(cells_.size());
				scalarField FNY(cells_.size());
				scalarField FNZ(cells_.size());
				scalarField FTX(cells_.size());
				scalarField FTY(cells_.size());
				scalarField FTZ(cells_.size());

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
				
				forAll(cells_, i)
				{
					if((CC[i].x()>=0.5) && (CC[i].x()<=1.5))

					{ 
					
					WX[i]= U[i].x();
					WY[i]= U[i].y() - OMEGAR;
					WZ[i]= 0;
					WMAG[i] = sqrt(WX[i]*WX[i] + WY[i]*WY[i] + WZ[i]*WZ[i]); 

					//Calculating WDOTN
					
					WDOTN[i] = WX[i]*NX + WY[i]*NY + WZ[i]*NZ;
					
////					//calculating WNX, WNY, WNZ, DEVLOC
////					
////					WNX[i] = WDOTN[i]*NX;
////					WNY[i] = WDOTN[i]*NY;
////					WNZ[i] = WDOTN[i]*NZ;
					DEVLOC[i] = asin(WDOTN[i]/WMAG[i]);
////					
////					//calculating WTX, WTY, WTZ, WTMAG
////					
////					
////					WTX[i] = WX[i] - WNX[i];
////					WTY[i] = WY[i] - WNY[i];
////					WTZ[i] = WZ[i] - WNZ[i];
////					WTMAG[i] = sqrt(WTX[i]*WTX[i] + WTY[i]*WTY[i] + WTZ[i]*WTZ[i]);
////
////					//calculating TX, TY, TZ
////					
////					TX[i] = WTX[i]/WTMAG[i];
////					TY[i] = WTY[i]/WTMAG[i];
////					TZ[i] = WTZ[i]/WTMAG[i];
					
					//calculating FNX, FNY, FNZ, FTX, FTY, FTZ
			
					FNX[i] = -PI*B*DEVLOC[i]*WMAG[i]*WMAG[i]*sin(DEVLOC[i]);
					FNY[i] = -PI*B*DEVLOC[i]*WMAG[i]*WMAG[i]*cos(DEVLOC[i]);
					FNZ[i] = -PI*B*DEVLOC[i]*WMAG[i]*WMAG[i]*NZ;
					//FNX[i] = -DEVLOC[i]*WMAG[i]*WMAG[i]*B*NX/R;
					//FNY[i] = -DEVLOC[i]*WMAG[i]*WMAG[i]*B*NY/R;
					//FNZ[i] = -DEVLOC[i]*WMAG[i]*WMAG[i]*B*NZ/R;

					//FTX[i] = DEVLOC[i]*WMAG[i]*WMAG[i]*B*TX[i]/R;
					//FTY[i] = DEVLOC[i]*WMAG[i]*WMAG[i]*B*TY[i]/R;
					//FTZ[i] = DEVLOC[i]*WMAG[i]*WMAG[i]*B*TZ[i]/R;

					//calculating MOMSRCX, MOMSRCY, MOMSRCZ
					
					MOMSRCX[i] = (FNX[i]);///cos(radTHETA));
					MOMSRCY[i] = (FNY[i]);///cos(radTHETA));
					MOMSRCZ[i] = (FNZ[i]);
					
					slope[i] = theta;
					angle[i] = - atan(U[i].y()/U[i].x());
					fangle[i] = (angle[i]*180/PI);
					deviation[i] = (DEVLOC[i]*180/PI);
					//
					//MOMSRCX[i] = (FTX[i] + FNX[i]);
					//MOMSRCY[i] = (FTY[i] + FNY[i]);
					//MOMSRCZ[i] = (FTZ[i] + FNZ[i]);
					//
					//adding source terms to the momentum equation
					 
					bodyForce[i] = vector(MOMSRCX[i], MOMSRCY[i], MOMSRCZ[i]);
					
					}
					
					else
					
					{
					
						MOMSRCX[i] = 0;
						MOMSRCY[i] = 0;
						MOMSRCZ[i] = 0;
					
						bodyForce[i] = vector(0, 0, 0);
					
					}
					//
					//std::cout << "x	y	z	slope	deviation	flow angle	" << std::endl;
					//
					
					//if (CC[i].y()=0.475)
					//{
					//if(mesh_.time().outputTime())
					//{
					//	std::cout << "	"	<< CC[i].x() << "	"<< CC[i].y() << "	"<< CC[i].z() <<
					//	"	" << theta <<
					//	"	" << DEVLOC[i] << 	
					//	"	" << std::endl;
					//	bodyForce.write();
					//}
					//}
		
					if(mesh_.time().outputTime())
					{
						bodyForce.write();
						// slope.write();
						// fangle.write();
						// deviation.write();
					}
					
			
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
    #line 37 "/home/palak/OpenFOAM/palak-6/run/Z1flatanalytical/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
//Known constant values: B, R, omega, theta
				const scalar B = 1;
				const scalar R = 1;
				const scalar omega = 1000;
				const scalar theta = 40;
				
				//Converting theta to radians
				scalar PI = constant::mathematical::pi;
				scalar radTHETA = (theta*PI/180.0);
				
				//Calculating OMEGAR, NX, NY, NZ
				scalar OMEGAR = omega*R;
				const scalar NX = sin(radTHETA);
				const scalar NY = cos(radTHETA);
				const scalar NZ = 0;
				
				//Initializing all fields
				
				scalarField x_coord(cells_.size());
				scalarField y_coord(cells_.size());
				
				scalarField WX(cells_.size());
				scalarField WY(cells_.size());
				scalarField WZ(cells_.size());
				scalarField WMAG(cells_.size());
				
				scalarField WDOTN(cells_.size());
				
				scalarField WNX(cells_.size());
				scalarField WNY(cells_.size());
				scalarField WNZ(cells_.size());
				scalarField DEVLOC(cells_.size());
				scalarField angle(cells_.size());
				
				scalarField WTX(cells_.size());
				scalarField WTY(cells_.size());
				scalarField WTZ(cells_.size());
				scalarField WTMAG(cells_.size());
				
				scalarField TX(cells_.size());
				scalarField TY(cells_.size());
				scalarField TZ(cells_.size());
				
				scalarField FNX(cells_.size());
				scalarField FNY(cells_.size());
				scalarField FNZ(cells_.size());
				scalarField FTX(cells_.size());
				scalarField FTY(cells_.size());
				scalarField FTZ(cells_.size());

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
				
				forAll(cells_, i)
				{
					if((CC[i].x()>=0.5) && (CC[i].x()<=1.5))

					{ 
					
					WX[i]= U[i].x();
					WY[i]= U[i].y() - OMEGAR;
					WZ[i]= 0;
					WMAG[i] = sqrt(WX[i]*WX[i] + WY[i]*WY[i] + WZ[i]*WZ[i]); 
					
					//Calculating WDOTN
					
					WDOTN[i] = WX[i]*NX + WY[i]*NY + WZ[i]*NZ;
					
////					//calculating WNX, WNY, WNZ, DEVLOC
////					
////					WNX[i] = WDOTN[i]*NX;
////					WNY[i] = WDOTN[i]*NY;
////					WNZ[i] = WDOTN[i]*NZ;
					DEVLOC[i] = asin(WDOTN[i]/WMAG[i]);
////					
////					//calculating WTX, WTY, WTZ, WTMAG
////					
////					
////					WTX[i] = WX[i] - WNX[i];
////					WTY[i] = WY[i] - WNY[i];
////					WTZ[i] = WZ[i] - WNZ[i];
////					WTMAG[i] = sqrt(WTX[i]*WTX[i] + WTY[i]*WTY[i] + WTZ[i]*WTZ[i]);
////
////					//calculating TX, TY, TZ
////					
////					TX[i] = WTX[i]/WTMAG[i];
////					TY[i] = WTY[i]/WTMAG[i];
////					TZ[i] = WTZ[i]/WTMAG[i];
					
					//calculating FNX, FNY, FNZ, FTX, FTY, FTZ
			
					FNX[i] = -PI*B*DEVLOC[i]*WMAG[i]*WMAG[i]*sin(DEVLOC[i]);
					FNY[i] = -PI*B*DEVLOC[i]*WMAG[i]*WMAG[i]*cos(DEVLOC[i]);
					FNZ[i] = -PI*B*DEVLOC[i]*WMAG[i]*WMAG[i]*NZ;
					//FNX[i] = -DEVLOC[i]*WMAG[i]*WMAG[i]*B*NX/R;
					//FNY[i] = -DEVLOC[i]*WMAG[i]*WMAG[i]*B*NY/R;
					//FNZ[i] = -DEVLOC[i]*WMAG[i]*WMAG[i]*B*NZ/R;

					//FTX[i] = DEVLOC[i]*WMAG[i]*WMAG[i]*B*TX[i]/R;
					//FTY[i] = DEVLOC[i]*WMAG[i]*WMAG[i]*B*TY[i]/R;
					//FTZ[i] = DEVLOC[i]*WMAG[i]*WMAG[i]*B*TZ[i]/R;

					//calculating MOMSRCX, MOMSRCY, MOMSRCZ
					
					MOMSRCX[i] = (FNX[i]);///cos(radTHETA));
					MOMSRCY[i] = (FNY[i]);///cos(radTHETA));
					MOMSRCZ[i] = (FNZ[i]);
					
					slope[i] = theta;
					angle[i] = - atan(U[i].y()/U[i].x());
					fangle[i] = (angle[i]*180/PI);
					deviation[i] = (DEVLOC[i]*180/PI);
					//
					//MOMSRCX[i] = (FTX[i] + FNX[i]);
					//MOMSRCY[i] = (FTY[i] + FNY[i]);
					//MOMSRCZ[i] = (FTZ[i] + FNZ[i]);
					//
					//adding source terms to the momentum equation
					 
					bodyForce[i] = vector(MOMSRCX[i], MOMSRCY[i], MOMSRCZ[i]);
					
					}
					
					else
					
					{
					
						MOMSRCX[i] = 0;
						MOMSRCY[i] = 0;
						MOMSRCZ[i] = 0;
					
						bodyForce[i] = vector(0, 0, 0);
					
					}
					//
					//std::cout << "x	y	z	slope	deviation	flow angle	" << std::endl;
					//
					
					//if (CC[i].y()=0.475)
					//{
					//if(mesh_.time().outputTime())
					//{
					//	std::cout << "	"	<< CC[i].x() << "	"<< CC[i].y() << "	"<< CC[i].z() <<
					//	"	" << theta <<
					//	"	" << DEVLOC[i] << 	
					//	"	" << std::endl;
					//	bodyForce.write();
					//}
					//}
		
					if(mesh_.time().outputTime())
					{
						bodyForce.write();
						// slope.write();
						// fangle.write();
						// deviation.write();
					}
					
			
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
    #line 276 "/home/palak/OpenFOAM/palak-6/run/Z1flatanalytical/constant/fvOptions.momentumSource.vectorCodedSourceCoeffs"
Pout<< "**codeSetValue**" << endl;
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

} // End namespace fv
// ************************************************************************* //


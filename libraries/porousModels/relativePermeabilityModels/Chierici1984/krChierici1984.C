/*---------------------------------------------------------------------------*\
  		  _______  ____    ____  ________  
 		 |_   __ \|_   \  /   _||_   __  | 
   		   | |__) | |   \/   |    | |_ \_| 
   		   |  ___/  | |\  /| |    |  _|    
    		  _| |_    _| |_\/_| |_  _| |_     
   		 |_____|  |_____||_____||_____|    
   	     Copyright (C) Toulouse INP, Pierre Horgue

License
    This file is part of porousMultiphaseFoam, an extension of OpenFOAM
    developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase 
    flows through porous media.

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

#include "krChierici1984.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{
defineTypeNameAndDebug(krChierici1984, 0);

addToRunTimeSelectionTable
(
    relativePermeabilityModel,
    krChierici1984,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePermeabilityModels::krChierici1984::krChierici1984
(
    const word& name,
    const dictionary& transportProperties,
    const volScalarField& Sb
)
    :
    relativePermeabilityModel(name, transportProperties,Sb),
    Smin_
    (
        IOobject(Sb_.name()+"min",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        transportProperties.lookupOrDefault(Sb_.name()+"min",dimensionedScalar(Sb_.name()+"min",dimless,0))
    ),
    Smax_
    (
        IOobject(Sb_.name()+"max",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        transportProperties.lookupOrDefault(Sb_.name()+"max",dimensionedScalar(Sb_.name()+"max",dimless,0))
    ),
    krChierici1984Coeffs_(transportProperties.subDict(typeName + "Coeffs")),
    cB_
    (
        IOobject("cb",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("cb",dimless,krChierici1984Coeffs_.lookupOrDefault<scalar>("cb",0))
    ),
    cM_
    (
        IOobject("cM",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("cM",dimless,krChierici1984Coeffs_.lookupOrDefault<scalar>("cM",0))
    ),
    cA_
    (
        IOobject("cA",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("cM",dimless,krChierici1984Coeffs_.lookupOrDefault<scalar>("cA",0))
    ),
    cL_
    (
        IOobject("cL",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("cL",dimless,krChierici1984Coeffs_.lookupOrDefault<scalar>("cL",0))
    ),
    kramax_
    (
        IOobject("kr"+Sb_.name()+"max",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("kr"+Sb_.name()+"max",dimless,krChierici1984Coeffs_.lookupOrDefault<scalar>("kr"+Sb_.name()+"max",1.0))
    ),
    krbmax_
    (
        IOobject
        (
            "kr"+Sb_.name()+"max",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("kr"+Sb_.name()+"max",dimless,krChierici1984Coeffs_.lookupOrDefault<scalar>("kr"+Sb_.name()+"max",1.0))
    )
{
    Se_ = (Sb_- Smin_)/(1 - Smax_-Smin_);
    if (gMin(cB_) <= 0 and gMin(cM_) <= 0 and gMin(cA_) <= 0 and gMin(cL_) <= 0)
    {
        FatalErrorIn
            (
                "in krChierici1984.C"
            )
            << "Relative permeability coefficient cB, cM, cA, cL equal or less than 0"
                << exit(FatalError);
    }
 
    Info << "Chierici1984 parameters for relative permeability model" << nl << "{" << endl;
    Info << "    cB ";
    if (cB_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(cB_).value() << endl;}
    Info << "    cM ";
    if (cA_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(cM_).value() << endl;}
    Info << "    cA ";
    if (cA_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(cA_).value() << endl;}
    Info << "    cL ";
    if (cL_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(cL_).value() << endl;}
    Info << "    Smax ";
    if (Smax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smax_).value() << endl;}
    Info << "    Smin ";
    if (Smin_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smin_).value() << endl;}
    Info << "    kramax ";
    if (kramax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(kramax_).value() << endl;}
    Info << "    krbmax ";
    if (krbmax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(krbmax_).value() << endl;}
    Info << "} \n" << endl;
}

// ************************************************************************* //

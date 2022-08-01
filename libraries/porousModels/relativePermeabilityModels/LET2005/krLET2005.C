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

#include "krLET2005.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{
defineTypeNameAndDebug(krLET2005, 0);

addToRunTimeSelectionTable
(
    relativePermeabilityModel,
    krLET2005,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePermeabilityModels::krLET2005::krLET2005
(
    const word& name,
    const dictionary& transportProperties,
    const volScalarField& Sb
)
    :
    relativePermeabilityModel(name, transportProperties,Sb),
    Smin_
    (
        IOobject
        (
            Sb_.name()+"min",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        transportProperties.lookupOrDefault(Sb_.name()+"min",dimensionedScalar(Sb_.name()+"min",dimless,0))
    ),
    Smax_
    (
        IOobject
        (
            Sb_.name()+"max",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        transportProperties.lookupOrDefault(Sb_.name()+"max",dimensionedScalar(Sb_.name()+"max",dimless,0))
    ),
    krLET2005Coeffs_(transportProperties.subDict(typeName + "Coeffs")),
    Lb_
    (
        IOobject("Lb",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("Lb",dimless,krLET2005Coeffs_.lookupOrDefault<scalar>("Lb",0))
    ),
    Eb_
    (
        IOobject("Eb",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("Eb",dimless,krLET2005Coeffs_.lookupOrDefault<scalar>("Eb",0))
    ),
    Tb_
    (
        IOobject("Tb",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("Tb",dimless,krLET2005Coeffs_.lookupOrDefault<scalar>("Tb",0))
    ),
    La_
    (
        IOobject("La",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("La",dimless,krLET2005Coeffs_.lookupOrDefault<scalar>("La",0))
    ),
    Ea_
    (
        IOobject("Ea",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("Ea",dimless,krLET2005Coeffs_.lookupOrDefault<scalar>("Ea",0))
    ),
    Ta_
    (
        IOobject("Ta",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("Ta",dimless,krLET2005Coeffs_.lookupOrDefault<scalar>("Ta",0))
    ),
    kramax_
    (
        IOobject("kr"+Sb_.name()+"max",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("kr"+Sb_.name()+"max",dimless,krLET2005Coeffs_.lookupOrDefault<scalar>("kr"+Sb_.name()+"max",1.0))
    ),
    krbmax_
    (
        IOobject
        (
            "kr"+Sb_.name()+"max",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("kr"+Sb_.name()+"max",dimless,krLET2005Coeffs_.lookupOrDefault<scalar>("kr"+Sb_.name()+"max",1.0))
    )
{
    Se_ = (Sb_- Smin_)/(1 - Smax_-Smin_);
    if (gMin(Lb_) <= 0 and gMin(Eb_) <= 0 and gMin(Tb_) <= 0 and gMin(La_) <= 0 and gMin(Ea_) <= 0 and gMin(Ta_) <= 0)
    {
        FatalErrorIn
            (
                "in krLET2005.C"
            )
            << "Relative permeability coefficients equal or less than 0"
                << exit(FatalError);
    }
 
    Info << "LET2005 parameters for relative permeability model" << nl << "{" << endl;
    Info << "    Lb ";
    if (Lb_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Lb_).value() << endl;}
    Info << "    Eb ";
    if (Eb_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Eb_).value() << endl;}
    Info << "    Tb ";
    if (Tb_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Tb_).value() << endl;}
    Info << "    La ";
    if (La_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(La_).value() << endl;}
    Info << "    Ea ";
    if (Ea_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Ea_).value() << endl;}
    Info << "    Ta ";
    if (Ta_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Ta_).value() << endl;}
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

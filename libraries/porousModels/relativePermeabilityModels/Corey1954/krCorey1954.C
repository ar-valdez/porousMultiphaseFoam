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

#include "krCorey1954.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{
defineTypeNameAndDebug(krCorey1954, 0);

addToRunTimeSelectionTable
(
    relativePermeabilityModel,
    krCorey1954,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePermeabilityModels::krCorey1954::krCorey1954
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
    krCorey1954Coeffs_(transportProperties.subDict(typeName + "Coeffs")),
    nb_
    (
        IOobject("nb",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("nb",dimless,krCorey1954Coeffs_.lookupOrDefault<scalar>("nb",0))
    ),
    na_
    (
        IOobject("na",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("na",dimless,krCorey1954Coeffs_.lookupOrDefault<scalar>("na",0))
    ),
    kramax_
    (
        IOobject("kr"+Sb_.name()+"max",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("kr"+Sb_.name()+"max",dimless,krCorey1954Coeffs_.lookupOrDefault<scalar>("kr"+Sb_.name()+"max",1.0))
    ),
    krbmax_
    (
        IOobject
        (
            "kr"+Sb_.name()+"max",Sb_.time().timeName(),Sb_.db(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),
        Sb.mesh(),
        dimensionedScalar("kr"+Sb_.name()+"max",dimless,krCorey1954Coeffs_.lookupOrDefault<scalar>("kr"+Sb_.name()+"max",1.0))
    )
{
    Se_ = (Sb_- Smin_)/(Smax_-Smin_);
    if (gMin(nb_) <= 0 and gMin(na_) <=0)
    {
        FatalErrorIn
            (
                "in krCorey1954.C"
            )
            << "Relative permeability coefficient n equal or less than 0"
                << exit(FatalError);
    }
 
    Info << "Corey1954 parameters for relative permeability model" << nl << "{" << endl;
    Info << "    nb ";
    if (nb_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(nb_).value() << endl;}
    Info << "    na ";
    if (na_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(na_).value() << endl;}
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

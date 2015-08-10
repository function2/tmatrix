// Copyright (C) 2012  Michael Seyfert <michael@codesand.org>
/*
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
Balance a chemical equation:
for example passing
Mg H2O -> Mg+2 OH-1 H2
on the command line will output the balanced equation (so the
atoms match, and the net charge matches.)

{compounds} -> {compounds}
compounds : {compound} [{compounds}]
compound : {element}[number][{element}...][+digit | -digit]

Whitespace between compounds is required.
*/

#include<iostream>
#include<map>
#include<vector>
#include<sstream>
#include<cctype>
#include<cerrno>
#include"las/matrix.h"
#include"las/mtypes.h"

using namespace std;
using namespace las;

// atom_t - These are in order by atomic number.
enum atom_t
{
	NOATOM,
	H,
	He,
	Li,
	Be,
	B,
	C,
	N,
	O,
	F,
	Ne,
	Na,
	Mg,
	Al,
	Si,
	P,
	S,
	Cl,
	Ar,
	K,
	Ca,
	Sc,
	Ti,
	V,
	Cr,
	Mn,
	Fe,
	Co,
	Ni,
	Cu,
	Zn,
	Ga,
	Ge,
	As,
	Se,
	Br,
	Kr,
	Rb,
	Sr,
	Y,
	Zr,
	Nb,
	Mo,
	Tc,
	Ru,
	Rh,
	Pd,
	Ag,
	Cd,
	In,
	Sn,
	Sb,
	Te,
	I,
	Xe,
	Cs,
	Ba,
	La,
	Ce,
	Pr,
	Nd,
	Pm,
	Sm,
	Eu,
	Gd,
	Tb,
	Dy,
	Ho,
	Er,
	Tm,
	Yb,
	Lu,
	Hf,
	Ta,
	W,
	Re,
	Os,
	Ir,
	Pt,
	Au,
	Hg,
	Tl,
	Pb,
	Bi,
	Po,
	At,
	Rn,
	Fr,
	Ra,
	Ac,
	Th,
	Pa,
	U,
	Np,
	Pu,
	Am,
	Cm,
	Bk,
	Cf,
	Es,
	Fm,
	Md,
	No,
	Lr,
	Rf,
	Db,
	Sg,
	Bh,
	Hs,
	Mt,
	Ds,
	Rg
};

const char* const ATOM_TO_STR[] =
{
	"NOATOM",
	"H",
	"He",
	"Li",
	"Be",
	"B",
	"C",
	"N",
	"O",
	"F",
	"Ne",
	"Na",
	"Mg",
	"Al",
	"Si",
	"P",
	"S",
	"Cl",
	"Ar",
	"K",
	"Ca",
	"Sc",
	"Ti",
	"V",
	"Cr",
	"Mn",
	"Fe",
	"Co",
	"Ni",
	"Cu",
	"Zn",
	"Ga",
	"Ge",
	"As",
	"Se",
	"Br",
	"Kr",
	"Rb",
	"Sr",
	"Y",
	"Zr",
	"Nb",
	"Mo",
	"Tc",
	"Ru",
	"Rh",
	"Pd",
	"Ag",
	"Cd",
	"In",
	"Sn",
	"Sb",
	"Te",
	"I",
	"Xe",
	"Cs",
	"Ba",
	"La",
	"Ce",
	"Pr",
	"Nd",
	"Pm",
	"Sm",
	"Eu",
	"Gd",
	"Tb",
	"Dy",
	"Ho",
	"Er",
	"Tm",
	"Yb",
	"Lu",
	"Hf",
	"Ta",
	"W",
	"Re",
	"Os",
	"Ir",
	"Pt",
	"Au",
	"Hg",
	"Tl",
	"Pb",
	"Bi",
	"Po",
	"At",
	"Rn",
	"Fr",
	"Ra",
	"Ac",
	"Th",
	"Pa",
	"U",
	"Np",
	"Pu",
	"Am",
	"Cm",
	"Bk",
	"Cf",
	"Es",
	"Fm",
	"Md",
	"No",
	"Lr",
	"Rf",
	"Db",
	"Sg",
	"Bh",
	"Hs",
	"Mt",
	"Ds",
	"Rg"
};

atom_t StrToAtom(string atom_str)
{
	static map<string,atom_t> am;
	if(am.empty()){
		// Initialize the atom map.
		am.insert(make_pair("H",H));
		am.insert(make_pair("He",He));
		am.insert(make_pair("Li",Li));
		am.insert(make_pair("Be",Be));
		am.insert(make_pair("B",B));
		am.insert(make_pair("C",C));
		am.insert(make_pair("N",N));
		am.insert(make_pair("O",O));
		am.insert(make_pair("F",F));
		am.insert(make_pair("Ne",Ne));
		am.insert(make_pair("Na",Na));
		am.insert(make_pair("Mg",Mg));
		am.insert(make_pair("Al",Al));
		am.insert(make_pair("Si",Si));
		am.insert(make_pair("P",P));
		am.insert(make_pair("S",S));
		am.insert(make_pair("Cl",Cl));
		am.insert(make_pair("Ar",Ar));
		am.insert(make_pair("K",K));
		am.insert(make_pair("Ca",Ca));
		am.insert(make_pair("Sc",Sc));
		am.insert(make_pair("Ti",Ti));
		am.insert(make_pair("V",V));
		am.insert(make_pair("Cr",Cr));
		am.insert(make_pair("Mn",Mn));
		am.insert(make_pair("Fe",Fe));
		am.insert(make_pair("Co",Co));
		am.insert(make_pair("Ni",Ni));
		am.insert(make_pair("Cu",Cu));
		am.insert(make_pair("Zn",Zn));
		am.insert(make_pair("Ga",Ga));
		am.insert(make_pair("Ge",Ge));
		am.insert(make_pair("As",As));
		am.insert(make_pair("Se",Se));
		am.insert(make_pair("Br",Br));
		am.insert(make_pair("Kr",Kr));
		am.insert(make_pair("Rb",Rb));
		am.insert(make_pair("Sr",Sr));
		am.insert(make_pair("Y",Y));
		am.insert(make_pair("Zr",Zr));
		am.insert(make_pair("Nb",Nb));
		am.insert(make_pair("Mo",Mo));
		am.insert(make_pair("Tc",Tc));
		am.insert(make_pair("Ru",Ru));
		am.insert(make_pair("Rh",Rh));
		am.insert(make_pair("Pd",Pd));
		am.insert(make_pair("Ag",Ag));
		am.insert(make_pair("Cd",Cd));
		am.insert(make_pair("In",In));
		am.insert(make_pair("Sn",Sn));
		am.insert(make_pair("Sb",Sb));
		am.insert(make_pair("Te",Te));
		am.insert(make_pair("I",I));
		am.insert(make_pair("Xe",Xe));
		am.insert(make_pair("Cs",Cs));
		am.insert(make_pair("Ba",Ba));
		am.insert(make_pair("La",La));
		am.insert(make_pair("Ce",Ce));
		am.insert(make_pair("Pr",Pr));
		am.insert(make_pair("Nd",Nd));
		am.insert(make_pair("Pm",Pm));
		am.insert(make_pair("Sm",Sm));
		am.insert(make_pair("Eu",Eu));
		am.insert(make_pair("Gd",Gd));
		am.insert(make_pair("Tb",Tb));
		am.insert(make_pair("Dy",Dy));
		am.insert(make_pair("Ho",Ho));
		am.insert(make_pair("Er",Er));
		am.insert(make_pair("Tm",Tm));
		am.insert(make_pair("Yb",Yb));
		am.insert(make_pair("Lu",Lu));
		am.insert(make_pair("Hf",Hf));
		am.insert(make_pair("Ta",Ta));
		am.insert(make_pair("W",W));
		am.insert(make_pair("Re",Re));
		am.insert(make_pair("Os",Os));
		am.insert(make_pair("Ir",Ir));
		am.insert(make_pair("Pt",Pt));
		am.insert(make_pair("Au",Au));
		am.insert(make_pair("Hg",Hg));
		am.insert(make_pair("Tl",Tl));
		am.insert(make_pair("Pb",Pb));
		am.insert(make_pair("Bi",Bi));
		am.insert(make_pair("Po",Po));
		am.insert(make_pair("At",At));
		am.insert(make_pair("Rn",Rn));
		am.insert(make_pair("Fr",Fr));
		am.insert(make_pair("Ra",Ra));
		am.insert(make_pair("Ac",Ac));
		am.insert(make_pair("Th",Th));
		am.insert(make_pair("Pa",Pa));
		am.insert(make_pair("U",U));
		am.insert(make_pair("Np",Np));
		am.insert(make_pair("Pu",Pu));
		am.insert(make_pair("Am",Am));
		am.insert(make_pair("Cm",Cm));
		am.insert(make_pair("Bk",Bk));
		am.insert(make_pair("Cf",Cf));
		am.insert(make_pair("Es",Es));
		am.insert(make_pair("Fm",Fm));
		am.insert(make_pair("Md",Md));
		am.insert(make_pair("No",No));
		am.insert(make_pair("Lr",Lr));
		am.insert(make_pair("Rf",Rf));
		am.insert(make_pair("Db",Db));
		am.insert(make_pair("Sg",Sg));
		am.insert(make_pair("Bh",Bh));
		am.insert(make_pair("Hs",Hs));
		am.insert(make_pair("Mt",Mt));
		am.insert(make_pair("Ds",Ds));
		am.insert(make_pair("Rg",Rg));
	}

	map<string,atom_t>::iterator I = am.find(atom_str);
	if(I == am.end())
		return NOATOM;
	return I->second;
}

//===========================================================
// Element
//
// Represents a chemical element.
//===========================================================
class Element
{
public:
	Element(atom_t a = NOATOM) : _atom(a){}
	Element(string atom_str)
		: _atom(StrToAtom(atom_str)){}

	inline string AtomName()const
	{return ATOM_TO_STR[_atom];}

	inline void SetAtom(atom_t a)
	{_atom = a;}

	inline bool None()const
	{return _atom == NOATOM;}

	bool operator<(const Element &e)const
	{
		return _atom < e._atom;
	}
	bool operator==(const Element &e)const
	{
		return _atom == e._atom;
	}

private:
	atom_t _atom;
};

ostream& operator<<(ostream &rs, const Element &e)
{
	return rs << e.AtomName();
}
//===========================================================

//===========================================================
// Compound
//
// Represents a chemical compound.
//===========================================================
class Compound
{
public:
	Compound(string cmpd_str = "");

	// Parse cmpd_str, there should be no whitespace in this string.
	bool SetCompound(string cmpd_str);

	// Get a string representation of the compound.
	string Str()const;

	// Return the compound.
	vector<pair<Element,int> > Get()const
	{return _cmpd;}

	inline int Charge()const
	{return _net_charge;}

	// Return the number of atoms in the compound.
	int NumAtoms(const Element &e)const;

	bool operator<(const Compound &c)const; // TODO
	bool operator==(const Compound &c)const; // TODO
private:
	// <Atom, # of atom in compound>
	vector<pair<Element,int> > _cmpd;

	// The net charge of the compound, by default this is zero.
	int _net_charge;
};

Compound::Compound(string cmpd_str)
 : _net_charge(0)
{
	SetCompound(cmpd_str);
}

bool Compound::SetCompound(string cmpd_str)
{
	// <atom, # in compound>
	map<int, int> atom_counter;

	bool reading_atom = false;
	int atom_start = 0;
	int idx = 0;
	for(;;)
	{
		if(reading_atom){
			// Go until we find a non lowercase letter.
			while(idx < (int)cmpd_str.size() && islower(cmpd_str[idx]))
				++idx;
			string atom_name = cmpd_str.substr(atom_start, idx - atom_start);
			//cerr << "Got " << atom_name << endl;
			reading_atom = false;
			// Check for a number after the atom.
			int num_start = idx;
			while(idx < (int)cmpd_str.size() && isdigit(cmpd_str[idx]))
				++idx;
			long num_atom = 1;
			char *endptr;
			errno = 0;
			if(idx != num_start){
				string tmp = cmpd_str.substr(num_start,idx - num_start);
				num_atom = strtol(tmp.c_str(),&endptr,0);
				if(*endptr != '\0' || errno != 0 || num_atom < 1) return false;
			}
			//cerr << "Got# " << num_atom << endl;

			// Add to counter
			atom_t atom = StrToAtom(atom_name);
			if(atom == NOATOM) return false; // Unknown atom: bad name.
			atom_counter[(int)atom] += num_atom;

			if(idx == (int)cmpd_str.size()) break;
		}else{
			// Go until we find an uppercase letter, which represents
			// the atom name.
			if(isupper(cmpd_str[idx])){
				reading_atom = true;
				atom_start = idx;
			}else if(cmpd_str[idx] == '+' || cmpd_str[idx] == '-'){
				break;
			}else{
				return false;
			}
			++idx;
		}
	}

	// All that remains is the charge at the end.
	_net_charge = 0;
	if(idx < (int)cmpd_str.size()){
		if(cmpd_str[idx] != '+' && cmpd_str[idx] != '-')
			return false;
		int num_start = idx++;
		while(idx < (int)cmpd_str.size() && isdigit(cmpd_str[idx]))
			++idx;
		long charge = 0;
		char *endptr;
		errno = 0;
		string tmp = cmpd_str.substr(num_start,idx - num_start);
		charge = strtol(tmp.c_str(),&endptr,0);
		if(*endptr != '\0' || errno != 0) return false;
		//cerr << "Got charge " << charge << endl;
		_net_charge = charge;
	}

	// Set _cmpd .
	_cmpd.clear();
	for(map<int,int>::const_iterator I = atom_counter.begin();I != atom_counter.end();++I)
	{
		_cmpd.push_back(make_pair(Element((atom_t)I->first),I->second));
	}
	// Make sure to sort it for comparisons later.
	sort(_cmpd.begin(), _cmpd.end());
	return true;
}

string Compound::Str()const
{
	ostringstream oss;
	for(size_t k = 0;k < _cmpd.size();++k){
		oss << _cmpd[k].first.AtomName();
		if(_cmpd[k].second != 1)
			oss << _cmpd[k].second;
	}
	// Add the charge.
	if(_net_charge > 0)
		oss << '+' << _net_charge;
	else if(_net_charge < 0)
		oss << _net_charge;
	return oss.str();
}

int Compound::NumAtoms(const Element &e) const
{
	int ret = 0;
	for(size_t k = 0;k < _cmpd.size();++k)
		if(_cmpd[k].first == e)
			ret += _cmpd[k].second;
	return ret;
}
//===========================================================

//===========================================================
// ChemEquation
//
// Represents a chemical equation.
//===========================================================
class ChemEquation
{
public:
	ChemEquation(string eq_str = "")
	{
		SetEq(eq_str);
	}

	bool SetEq(string eq_str);

	inline bool Invalid()const
	{return _sides[0].empty() || _sides[1].empty();}

	string Str()const;

private:
	typedef vector<pair<Compound,bigq> > side_t;
	side_t _sides[2]; // left / right side of the equation.

	// Parse compounds and put them in lst.
	// Compounds are separated by whitespace.
	bool ParseCompounds(side_t &lst, string cmpds);

	// Balance the equation.
	bool Balance();
};

bool ChemEquation::SetEq(string eq_str)
{
	// Find the center of the equation.
	// This splits the left from the right side.
	int center = eq_str.find("->");
	if(center == -1) return false;

	string left = eq_str.substr(0,center);
	string right= eq_str.substr(center+2);

	if(!ParseCompounds(_sides[0],left))
		return false;
	if(!ParseCompounds(_sides[1],right))
		return false;

	if(!Balance()){
		_sides[0].clear();
		_sides[1].clear();
		return false;
	}
	return true;
}

bool ChemEquation::ParseCompounds(side_t &lst, string cmpds)
{
	istringstream iss(cmpds);
	string t;
	while(iss >> t){
		Compound c;
		if(!c.SetCompound(t))
			return false;
		lst.push_back(make_pair(c,bigq(0)));
	}
	return true;
}

bool ChemEquation::Balance()
{
	if(Invalid()){
		cerr << "Invalid equation, unable to balance.\n";
		return false;
	}

	map<Element,int> atoms[2]; // Atoms on each side of the equation.

	for(size_t sidx = 0;sidx < 2;++sidx)
		for(size_t idx = 0;idx < _sides[sidx].size();++idx){
			vector<pair<Element,int> > p = _sides[sidx][idx].first.Get();
			for(size_t k = 0;k < p.size();++k){
				atoms[sidx][p[k].first] += p[k].second;
			}
		}

	// Make sure the atoms match, otherwise it is impossible to balance.
	for(map<Element,int>::const_iterator I = atoms[0].begin(); I != atoms[0].end();++I){
		if(atoms[1].find(I->first) == atoms[1].end()){
			cerr << "Atom " << I->first << " does not exist on other side. Unable to balance the equation.\n";
			return false;
		}
	}
	for(map<Element,int>::const_iterator I = atoms[1].begin(); I != atoms[1].end();++I){
		if(atoms[0].find(I->first) == atoms[0].end()){
			cerr << "Atom " << I->first << " does not exist on other side. Unable to balance the equation.\n";
			return false;
		}
	}

	/*
	  Each row represents an atom (element) in the equation.
	  The number of columns is equal to the number of separate compounds/elements
	   in the equation.
	  Each row will be filled with a relation between the left and right side.

	  Example:
	  Mg H2O -> Mg+2 OH-1 H2
	  A  B      C    D    E

	  Mg [ -1  0  1  0  0 ] * [A B C D E]^T = 0
	   H [  0 -2  0  1  2 ]
		O [  0 -1  0  1  0 ]
	  +- [  0  0  2 -1  0 ] // The charges must match too.
	 */

	// Number of rows is the number of atoms, plus one for the charge relation.
	// Number of cols is the number of unknowns (compound count).
	mxqr A(atoms[0].size() + 1, _sides[0].size() + _sides[1].size());

	// For each element.
	int row = 0;
	for(map<Element,int>::const_iterator I = atoms[0].begin();I != atoms[0].end();++I,++row){
		Element e = I->first;
		int col = 0;
		for(size_t sidx = 0;sidx < 2;++sidx)
			for(size_t k = 0;k < _sides[sidx].size();++k,++col)
				A(row,col) = _sides[sidx][k].first.NumAtoms(e) * (sidx ? 1 : -1);
	}

	// Add the charge numbers relation.
	int col = 0;
	for(size_t sidx = 0;sidx < 2;++sidx)
		for(size_t k = 0;k < _sides[sidx].size();++k,++col)
			A(row,col) = _sides[sidx][k].first.Charge() * (sidx ? 1 : -1);

	cerr << "=DEBUG= EQ is:\n";
	MatrixOutput(cerr,A);
	cerr << endl;

	// Find the nullspace of the equation.
	mxqr R; // Reduced row echelon form.
	RRef(A,R);
	mxqr N;
	Nullspace(R,N);

	cerr << "=DEBUG= Nullspace is:\n";
	MatrixOutput(cerr,N);
	cerr << endl;

	if(N.size2() == 0){
		cerr << "No nullspace? Unable to balance.";
		return false;
	}

	// Get the LCM of the denominators, so that our final result is all integers.
	bigint lcm_val = 1;
	col = 0;
	for(row = 0;row < (int)N.size1();++row){
		bigint tmp = N(row,col).get_den();
		bigint old_lcm = lcm_val;
		mpz_lcm(lcm_val.get_mpz_t(),old_lcm.get_mpz_t(),tmp.get_mpz_t());
	}

	row = 0;
	for(size_t sidx = 0;sidx < 2;++sidx)
		for(size_t k = 0;k < _sides[sidx].size();++k,++row)
			_sides[sidx][k].second = N(row,col) * lcm_val;

	return true;
}

string ChemEquation::Str()const
{
	ostringstream oss;

	bool first = true;
	for(size_t sidx = 0;sidx < 2;++sidx){
		for(size_t k = 0;k < _sides[sidx].size();++k){
			if(!first) oss << ' ';
			else first = false;
			bigq num = _sides[sidx][k].second;
			//if(num != 1)
			oss << _sides[sidx][k].second;
			//oss << _sides[sidx][k].first.Str();
		}
		if(sidx == 0)
			;//oss << " -> ";
	}
	return oss.str();
}
//===========================================================

int main(int argc, char**argv)
{
	ChemEquation eq(argv[1]);
	if(!eq.Invalid())
		cout << eq.Str() << endl;
}

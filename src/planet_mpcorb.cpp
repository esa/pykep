#include<fstream>
#include<boost/algorithm/string.hpp>

#include"planet_mpcorb.h"
#include"exceptions.h"

static const int mpcorb_format[8][2] =
{
	{92,11},	// a (AU)
	{70,9},		// e
	{59,9},		// i (deg)
	{48,9},		// Omega (deg)
	{37,9},		// omega (deg)
	{26,9},		// M (deg)
	{20,5},		// Epoch (shitty format)
	{166,28}	// Asteroid readable name
};

namespace kep_toolbox{

planet_mpcorb::planet_mpcorb(const std::string& line)
{
	std::string linecopy(line);
	boost::algorithm::to_lower(linecopy);
	array6D elem;
	std::string tmp;
	//read keplerian elements from MPCORB.DAT
	for (int i = 0; i < 6; ++i) {
		tmp.clear();
		tmp.append(&linecopy[mpcorb_format[i][0]],mpcorb_format[i][1]);
		boost::algorithm::trim(tmp);
		elem[i] = boost::lexical_cast<double>(tmp);
	}
	// Converting orbital elements to the dictatorial PaGMO units.
	elem[0] *= ASTRO_AU;
	for (int i = 2; i < 6; ++i) {
		elem[i] *= ASTRO_DEG2RAD;
	}
	// Deal with MPCORB data format
	tmp.clear();
	tmp.append(&linecopy[mpcorb_format[6][0]],mpcorb_format[6][1]);
	boost::algorithm::trim(tmp);
	kep_toolbox::epoch epoch(packed_date2epoch(tmp));

	// Record asteroid name.
	tmp.clear();
	tmp.append(&linecopy[mpcorb_format[7][0]],mpcorb_format[7][1]);
	boost::algorithm::trim(tmp);
	build_planet(epoch,elem,ASTRO_MU_SUN,100,100,100,tmp);
}


epoch planet_mpcorb::packed_date2epoch(std::string in) {
	if (in.size()!=5) {
		throw_value_error("mpcorb data format requires 5 characters.");
	}
	boost::algorithm::to_lower(in);
	boost::gregorian::greg_year anno = packed_date2number(in[0]) * 100 + boost::lexical_cast<int>(std::string(&in[1],&in[3]));
	boost::gregorian::greg_month mese = packed_date2number(in[3]);
	boost::gregorian::greg_day giorno = packed_date2number(in[4]);
	return epoch(anno,mese,giorno);
}
// Convert mpcorb packed dates convention into number. (lower case assumed)
// TODO: check locale ASCII.
int planet_mpcorb::packed_date2number(char c)
{
	return static_cast<int>(c) - (boost::algorithm::is_alpha()(c) ? 87 : 48);
}

planet_ptr planet_mpcorb::clone() const
{
	return planet_ptr(new planet_mpcorb(*this));
}

} //namespace

// Serialization code
BOOST_CLASS_EXPORT_IMPLEMENT(kep_toolbox::planet_mpcorb);
// Serialization code (END)

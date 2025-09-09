/**
 * Creates needed tables
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date nov-2015
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */


#include <regex>


// ============================================================================


t_cplx get_number(std::string str)
{
	tl::string_rm<std::string>(str, "(", ")");

	if(tl::string_rm<std::string>(str, "<i>", "</i>"))
	{	// complex number
		std::size_t iSign = str.find_last_of("+-");
		str.insert(iSign, ", ");
		str.insert(0, "(");
		str.push_back(')');
	}

	tl::trim(str);
	if(str == "---") str = "";


	t_cplx c;
	std::istringstream istr(str);
	istr >> c;

	t_real dR = c.real(), dI = c.imag();
	tl::set_eps_0(dR); tl::set_eps_0(dI);
	c.real(dR); c.imag(dI);

	return c;
}


bool get_abundance_or_hl(const std::string& _str, t_real& dAbOrHL)
{
	bool bIsHL = (_str.find("a") != std::string::npos);
	std::string str = tl::remove_chars(_str, std::string("()a"));

	dAbOrHL = get_number(str).real();
	if(!bIsHL) dAbOrHL /= t_real(100.);
	return !bIsHL;
}


bool gen_scatlens()
{
	std::ifstream ifstr("tmp/scatlens.html");
	if(!ifstr)
	{
		tl::log_err("Cannot open \"tmp/scatlens.html\".");
		return false;
	}

	bool bTableStarted = 0;
	std::string strTable;
	while(!ifstr.eof())
	{
		std::string strLine;
		std::getline(ifstr, strLine);

		if(!bTableStarted)
		{
			std::vector<std::string> vecHdr;
			tl::get_tokens_seq<std::string, std::string, std::vector>(strLine, "<th>", vecHdr, 0);
			if(vecHdr.size() < 9)
				continue;
			bTableStarted = 1;
		}
		else
		{
			// at end of table?
			if(tl::str_to_lower(strLine).find("/table") != std::string::npos)
				break;

			strTable += strLine;
		}
	}
	ifstr.close();



	std::vector<std::string> vecRows;
	tl::get_tokens_seq<std::string, std::string, std::vector>(strTable, "<tr>", vecRows, 0);


	tl::Prop<std::string> prop;
	prop.SetSeparator('.');
	prop.Add("scatlens.source", "Scattering lengths and cross-sections extracted from NIST table"
		" (which itself is based on <a href=\"http://dx.doi.org/10.1080/10448639208218770\">this paper</a>).");
	prop.Add("scatlens.source_url", "https://www.ncnr.nist.gov/resources/n-lengths/list.html");

	unsigned int iAtom = 0;
	for(const std::string& strRow : vecRows)
	{
		if(strRow.length() == 0)
			continue;

		std::ostringstream ostr;
		ostr << "scatlens.atom_" << iAtom;
		std::string strAtom = ostr.str();


		std::vector<std::string> vecCol;
		tl::get_tokens_seq<std::string, std::string, std::vector>(strRow, "<td>", vecCol, 0);
		if(vecCol.size() < 9)
		{
			tl::log_warn("Invalid number of table entries in row \"", strRow, "\".");
			continue;
		}

		std::string strName = vecCol[1];
		tl::trim(strName);
		if(strName == "") continue;

		t_cplx cCoh = get_number(vecCol[3]);
		t_cplx cIncoh = get_number(vecCol[4]);
		t_real dXsecCoh = get_number(vecCol[5]).real();
		t_real dXsecIncoh = get_number(vecCol[6]).real();
		t_real dXsecScat = get_number(vecCol[7]).real();
		t_real dXsecAbsTherm = get_number(vecCol[8]).real();

		t_real dAbOrHL = t_real(0);
		bool bAb = get_abundance_or_hl(vecCol[2], dAbOrHL);

		prop.Add(strAtom + ".name", strName);
		prop.Add(strAtom + ".coh", tl::var_to_str(cCoh, g_iPrec));
		prop.Add(strAtom + ".incoh", tl::var_to_str(cIncoh, g_iPrec));

		prop.Add(strAtom + ".xsec_coh", tl::var_to_str(dXsecCoh, g_iPrec));
		prop.Add(strAtom + ".xsec_incoh", tl::var_to_str(dXsecIncoh, g_iPrec));
		prop.Add(strAtom + ".xsec_scat", tl::var_to_str(dXsecScat, g_iPrec));
		prop.Add(strAtom + ".xsec_abs", tl::var_to_str(dXsecAbsTherm, g_iPrec));

		if(bAb)
			prop.Add(strAtom + ".abund", tl::var_to_str(dAbOrHL, g_iPrec));
		else
			prop.Add(strAtom + ".hl", tl::var_to_str(dAbOrHL, g_iPrec));

		++iAtom;
	}

	prop.Add("scatlens.num_atoms", tl::var_to_str(iAtom));


	if(!prop.Save("res/data/scatlens.xml.gz"))
	{
		tl::log_err("Cannot write \"res/data/scatlens.xml.gz\".");
		return false;
	}
	return true;
}



// ============================================================================


bool gen_magformfacts()
{
	tl::Prop<std::string> propOut;
	propOut.SetSeparator('.');
	propOut.Add("magffacts.source", "Magnetic form factor coefficients extracted from ILL table.");
	propOut.Add("magffacts.source_url", "https://www.ill.eu/sites/ccsl/ffacts/");

	std::size_t iAtom=0;
	std::set<std::string> setAtoms;

	std::vector<std::string> vecFiles =
		{"tmp/j0_1.html", "tmp/j0_2.html",
		"tmp/j0_3.html" , "tmp/j0_4.html",

		"tmp/j2_1.html", "tmp/j2_2.html",
		"tmp/j2_3.html", "tmp/j2_4.html",};

	for(std::size_t iFile=0; iFile<vecFiles.size(); ++iFile)
	{
		const std::string& strFile = vecFiles[iFile];
		std::string strJ = iFile < 4 ? "j0" : "j2";

		// switching to j2 files
		if(iFile==4)
		{
			iAtom = 0;
			setAtoms.clear();
		}

		std::ifstream ifstr(strFile);
		if(!ifstr)
		{
			tl::log_err("Cannot open \"", strFile, "\".");
			return false;
		}

		std::string strTable;
		bool bTableStarted=0;
		while(!ifstr.eof())
		{
			std::string strLine;
			std::getline(ifstr, strLine);
			std::string strLineLower = tl::str_to_lower(strLine);

			if(bTableStarted)
			{
				strTable += strLine + "\n";
				if(strLineLower.find("</table") != std::string::npos)
				break;
			}
			else
			{
				std::size_t iPos = strLineLower.find("<table");
				if(iPos != std::string::npos)
				{
					std::string strSub = strLine.substr(iPos);
					strTable += strSub + "\n";

					bTableStarted = 1;
				}
			}
		}
		if(strTable.length() == 0)
		{
			tl::log_err("Invalid table: \"", strFile, "\".");
			return 0;
		}

		// removing attributes
		std::basic_regex<char> rex("<([A-Za-z]*)[A-Za-z0-9\\=\\\"\\ ]*>", std::regex::ECMAScript);
		strTable = std::regex_replace(strTable, rex, "<$1>");

		tl::find_all_and_replace<std::string>(strTable, "<P>", "");
		tl::find_all_and_replace<std::string>(strTable, "<p>", "");


		std::istringstream istrTab(strTable);
		tl::Prop<std::string> prop;
		prop.Load(istrTab, tl::PropType::XML);
		const auto& tab = prop.GetProp().begin()->second;

		auto iter = tab.begin(); ++iter;
		for(; iter!=tab.end(); ++iter)
		{
			auto iterElem = iter->second.begin();
			std::string strElem = iterElem++->second.data();
			tl::trim(strElem);
			if(*strElem.rbegin() == '0')
				strElem.resize(strElem.length()-1);
			else
				strElem += "+";

			if(setAtoms.find(strElem) != setAtoms.end())
			{
				tl::log_warn("Atom ", strElem, " already in set. Ignoring.");
				continue;
			}
			setAtoms.insert(strElem);

			t_real dA = tl::str_to_var<t_real>(iterElem++->second.data());
			t_real da = tl::str_to_var<t_real>(iterElem++->second.data());
			t_real dB = tl::str_to_var<t_real>(iterElem++->second.data());
			t_real db = tl::str_to_var<t_real>(iterElem++->second.data());
			t_real dC = tl::str_to_var<t_real>(iterElem++->second.data());
			t_real dc = tl::str_to_var<t_real>(iterElem++->second.data());
			t_real dD = tl::str_to_var<t_real>(iterElem->second.data());

			std::string strA = tl::var_to_str(dA, g_iPrec) + "; " +
				tl::var_to_str(dB, g_iPrec) + "; " +
				tl::var_to_str(dC, g_iPrec) + "; " +
				tl::var_to_str(dD, g_iPrec);
			std::string stra = tl::var_to_str(da, g_iPrec) + "; " +
				tl::var_to_str(db, g_iPrec) + "; " +
				tl::var_to_str(dc, g_iPrec);

			std::ostringstream ostrAtom;
			ostrAtom << "magffacts." + strJ + ".atom_" << iAtom;
			propOut.Add(ostrAtom.str() + ".name", strElem);
			propOut.Add(ostrAtom.str() + ".A", strA);
			propOut.Add(ostrAtom.str() + ".a", stra);
			++iAtom;
		}
	}

	propOut.Add("magffacts.num_atoms", tl::var_to_str(iAtom));

	if(!propOut.Save("res/data/magffacts.xml.gz"))
	{
		tl::log_err("Cannot write \"res/data/magffacts.xml.gz\".");
		return false;
	}

	return true;
}



// ============================================================================

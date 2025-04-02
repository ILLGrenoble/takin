/**
 * random numbers
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 16-aug-2013
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#include "rand.h"
#include "../log/log.h"

#include <cstdlib>
#include <exception>
#include <time.h>
#include <sys/time.h>


namespace tl {

// thread-global engine to use if a thread-local has not been seeded
static std::mt19937 g_randeng_fallback;
static bool g_bFallbackInited = false;

static thread_local std::mt19937/*_64*/ g_randeng;
static thread_local bool g_bHasEntropy = false;
static thread_local bool g_bIsSeeded = false;


std::mt19937& get_randeng()
{
	if(g_bIsSeeded)
		return g_randeng;
	return g_randeng_fallback;
}


unsigned int get_rand_seed()
{
	// seed 0: random device
	unsigned int uiSeed0 = 0;
	try
	{
		std::random_device rnd;
		g_bHasEntropy = (rnd.entropy() != 0);
		uiSeed0 = rnd();
	}
	catch(const std::exception& ex)
	{
		log_debug(ex.what());
		uiSeed0 = 0;
	}

	// seed 1: time based
	struct timeval timev;
	gettimeofday(&timev, 0);
	unsigned int uiSeed1 = timev.tv_sec ^ timev.tv_usec;

	// total seed
	unsigned int uiSeed = uiSeed0 ^ uiSeed1;
	return uiSeed;
}


void init_rand()
{
	init_rand_seed(get_rand_seed());
}


void init_rand_seed(unsigned int uiSeed, bool debug_msg)
{
	std::string strEntr;

	if(debug_msg)
	{
		if(!g_bHasEntropy)
			strEntr = ", but entropy is zero";
		log_debug("Random seed: ", uiSeed, strEntr, ".");
	}

	srand(uiSeed);
	g_randeng = std::mt19937/*_64*/(uiSeed);
	g_bIsSeeded = true;

	// copy first engine to thread-global fallback
	if(!g_bFallbackInited)
	{
		g_randeng_fallback = g_randeng;
		g_bFallbackInited = true;
	}
}


unsigned int simple_rand(unsigned int iMax)
{
	return rand() % iMax;
}


template<>
bool rand_minmax(bool tMin, bool tMax)
{
	using t_rand = tl::_rand_int<unsigned char>;

	t_rand rand;
	return rand(tMin, tMax) == 1;
}

}

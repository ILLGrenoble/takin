/**
 * GL drawing
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 22-dec-2014
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

/**
 * Freetype rendering under OpenGL, inspired by:
 * https://en.wikibooks.org/wiki/OpenGL_Programming/Modern_OpenGL_Tutorial_Text_Rendering_01
 * https://en.wikibooks.org/wiki/OpenGL_Programming/Modern_OpenGL_Tutorial_Text_Rendering_02
 */

#ifndef __GL_FONT_H__
#define __GL_FONT_H__

#include "gl.h"

#define DEF_FONT "/usr/share/fonts/dejavu/DejaVuSansMono.ttf"
#define DEF_FONT_SIZE 12

#include <ft2build.h>
#include FT_FREETYPE_H

namespace tl {
class FontMap
{
	protected:
		bool m_bOk = false;
		FT_Library m_ftLib = 0;
		FT_Face m_ftFace = 0;

		int m_iCharsPerLine = 0, m_iLines = 0;
		int m_iTileH = 0, m_iTileW = 0;
		int m_iPadH = DEF_FONT_SIZE/4, m_iPadW = DEF_FONT_SIZE/4;
		int m_iLargeW = 0, m_iLargeH = 0;
		unsigned char *m_pcLarge = nullptr;

		static const std::string m_strCharset;
		using t_offsmap = std::unordered_map<
			typename std::string::value_type,
			std::pair<int, int>>;
		t_offsmap m_mapOffs;

	protected:
		void UnloadFont();

	static void draw_tile(unsigned char* pcBuf,
		unsigned int iW, unsigned int iH,
		unsigned int iTileW, unsigned int iTileH,
		unsigned int iPosX, unsigned int iPosY,
		const unsigned char* pcGlyph,
		unsigned int iGlyphXOffs, unsigned int iGlyphYOffs,
		unsigned int iGlyphW, unsigned int iGlyphH);

	public:
		FontMap(const char* pcFont, int iSize);
		FontMap();
		virtual ~FontMap();

		bool LoadFont(const char* pcFont, int iSize);
		bool LoadFont(FT_Face ftFace, int iSize);

		void dump(std::ostream& ostr) const;
		void dump(std::ostream& ostr, const std::pair<int,int>& pair) const;
		virtual bool IsOk() const { return m_bOk; }

		const unsigned char* GetBuffer() const { return m_pcLarge; }
		std::pair<int, int> GetOffset(char ch) const;
};


template<class T = GLdouble>
class GlFontMap : public FontMap
{
	protected:
		bool m_bOk = false;
		GLuint m_tex = 0;
		T m_dScale = 0.01;

	protected:
		bool CreateFontTexture();

	public:
		GlFontMap() = delete;
		GlFontMap(const char* pcFont = DEF_FONT, int iSize = DEF_FONT_SIZE);
		GlFontMap(FT_Face ftFace, int iSize=DEF_FONT_SIZE);
		virtual ~GlFontMap();

		virtual bool IsOk() const override { return m_bOk && FontMap::IsOk(); }

		void BindTexture();
		void DrawText(T dX, T dY, const std::string& str, bool bCenter = true);
		void DrawText(T dX, T dY, T dZ, const std::string& str, bool bCenter = true);

		void SetScale(T dScale) { m_dScale = dScale; }
};

// --------------------------------------------------------------------------------

}

#endif

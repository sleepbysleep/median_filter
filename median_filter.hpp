/*
  Copyright (C) 2014 Hoyoung Lee

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

#pragma once

#include <cstdio>
#include <cstdlib>
//#include <fstream>
#include <iostream>
#include <cassert>
//#include <string>
//#include <array>
//#include <vector>
//#include <iterator>

#include <cpixmap.hpp>
#include <cchunk.hpp>

template <typename T>
// ascending or descending sorting by partial selection
inline T getKthElement(bool ascending, T *sorting_bin, size_t n, size_t k)
{
  for (size_t i = 0; i <= k; i++) {
    size_t index = i;
    T value = sorting_bin[i];    
    
    for (size_t j = i+1; j < n; j++) {
      if (ascending) { // ascending order
	if (sorting_bin[j] < value) {
	  index = j;
	  value = sorting_bin[j];
	}
      } else { // descending order
	if (sorting_bin[j] > value) {
	  index = j;
	  value = sorting_bin[j];
	}
      }
    }
    
    if (index != i) {
      value = sorting_bin[i];
      sorting_bin[i] = sorting_bin[index];
      sorting_bin[index] = value;
    }
  }
}

template <typename T>
void filterMedianKernel(cpixmap<T>& dst, const cpixmap<T>& src, const cregion<size_t>& window)
{
  assert(dst.isMatched(src));

  int loffset = window.getLeftHalf();
  int roffset = window.getRightHalf();
  int hpadding = window.getWidth()>>1;
  
  int uoffset = window.getUpHalf();
  int doffset = window.getDownHalf();
  int vpadding = window.getHeight()>>1;
  
  int n = window.getWidth() * window.getHeight();
  T *sorting_bin = new T[n];

  cslice<T> linebuf(src, 1/*lines*/, hpadding, vpadding);
  for (size_t z = 0; z < src.getBands(); z++) {
    linebuf.draftSlice(src, z);
    for (size_t y = 0; y < src.getBands(); y++) {
      T *dstline = dst.getLine(y, z);      
#pragma omp parallel for
      for (size_t x = 0; x < src.getWidth(); x++) {
	//std::array<T, n> sorting_bin;
	int k = 0;
	for (int i = y+uoffset; i < y+doffset; i++) {
	  for (int j = x+loffset; j < x+roffset; j++) {
	    sorting_bin[k++] = linebuf(i, j);
	  }
	}
	//std::nth_element(sorting_bin.begin(), sorting_bin.begin()+sorting_bin.size()/2, sorting_bin.end());
	//*(dstline+x) = sorting_bin[sorting_bin.size()/2];
	*(dstline+x) = getKthElement(true, sorting_bin, k, k/2);

      }
      linebuf.shiftSlice(1, src, z);      
    }
  }
  delete [] sorting_bin;
}

template <typename T>
void filterMedian3x3Kernel(cpixmap<T>& dst, const cpixmap<T>& src)
{
  assert(dst.isMatched(src));
  
  window3x3_frame<T> linebuf(src);
  for (size_t z = 0; z < src.getBands(); z++) {
    linebuf.draftFrame(src, z);
    for (size_t y = 0; y < src.getHeight(); y++) {
      T *dstline = dst.getLine(y, z);
#pragma omp parallel for
      for (size_t x = 0; x < src.getWidth(); x++) {
	T sorting_bin[9];	
	sorting_bin[0] = linebuf(y-1, x-1);
	sorting_bin[1] = linebuf(y-1, x);
	sorting_bin[2] = linebuf(y-1, x+1);
	sorting_bin[3] = linebuf(y, x-1);
	sorting_bin[4] = linebuf(y, x);
	sorting_bin[5] = linebuf(y, x+1);
	sorting_bin[6] = linebuf(y+1, x-1);
	sorting_bin[7] = linebuf(y+1, x);
	sorting_bin[8] = linebuf(y+1, x+1);
	for (size_t i = 0; i <= 4; i++) {
	  size_t index = i;
	  T value = sorting_bin[i];
	  for (size_t j = i+1; j < 9; j++) {
	    if (sorting_bin[j] < value) {
	      value = sorting_bin[j];
	      index = j;
	    }
	  }
	  if (index != i) {
	    value = sorting_bin[i];
	    sorting_bin[i] = sorting_bin[index];
	    sorting_bin[index] = value;	    
	  }
	}
	*(dstline+x) = sorting_bin[4];
      }
      linebuf.shiftFrame(src, z);      
    }
  }
}

template <typename T>
void filterMedian5x5Kernel(cpixmap<T>& dst, const cpixmap<T>& src)
{
  assert(dst.isMatched(src));

  window5x5_frame<T> linebuf(src);
  for (size_t z = 0; z < src.getBands(); z++) {
    linebuf.draftFrame(src, z);
    for (size_t y = 0; y < src.getHeight(); y++) {
      T *dstline = dst.getLine(y, z);
#pragma omp parallel for
      for (size_t x = 0; x < src.getWidth(); x++) {
	T sorting_bin[25];	
	sorting_bin[0] = linebuf(y-2, x-2);
	sorting_bin[1] = linebuf(y-2, x-1);
	sorting_bin[2] = linebuf(y-2, x);
	sorting_bin[3] = linebuf(y-2, x-1);
	sorting_bin[4] = linebuf(y-2, x+2);
	sorting_bin[5] = linebuf(y-1, x-2);
	sorting_bin[6] = linebuf(y-1, x-1);
	sorting_bin[7] = linebuf(y-1, x);
	sorting_bin[8] = linebuf(y-1, x+1);
	sorting_bin[9] = linebuf(y-1, x+2);
	sorting_bin[10] = linebuf(y, x-2);
	sorting_bin[11] = linebuf(y, x-1);
	sorting_bin[12] = linebuf(y, x);
	sorting_bin[13] = linebuf(y, x+1);
	sorting_bin[14] = linebuf(y, x+2);
	sorting_bin[15] = linebuf(y+1, x-2);
	sorting_bin[16] = linebuf(y+1, x-1);
	sorting_bin[17] = linebuf(y+1, x);
	sorting_bin[18] = linebuf(y+1, x+1);
	sorting_bin[19] = linebuf(y+1, x+2);
	sorting_bin[20] = linebuf(y+2, x-2);
	sorting_bin[21] = linebuf(y+2, x-1);
	sorting_bin[22] = linebuf(y+2, x);
	sorting_bin[23] = linebuf(y+2, x+1);
	sorting_bin[24] = linebuf(y+2, x+2);

	for (size_t i = 0; i <= 12; i++) {
	  size_t index = i;
	  T value = sorting_bin[i];
	  for (size_t j = i+1; j < 25; j++) {
	    if (sorting_bin[j] < value) {
	      value = sorting_bin[j];
	      index = j;
	    }
	  }
	  if (index != i) {
	    value = sorting_bin[i];
	    sorting_bin[i] = sorting_bin[index];
	    sorting_bin[index] = value;	    
	  }
	}
	*(dstline+x) = sorting_bin[12];
      }
      linebuf.shiftFrame(src, z);      
    }
  }
}

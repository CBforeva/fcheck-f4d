#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>


// Version identifier: year, month, day, release number
const int   VERSION_ID = 0x17031500;
//const int VERSION_ID = 0x17022700;
//const int VERSION_ID = 0x17013100;
//const int VERSION_ID = 0x16120200;

// FROM FELIXDEFS:
#define BLOCK_BYTES                1024
#define BLOCK_SHORTS              (BLOCK_BYTES/2)
#define BLOCK_HEADER_BYTES         4

#define BLOCK_ID                   0xABCD

#define BLOCK_SEQNR_MASK           0xF800
#define BLOCK_SEQNR_SHIFT          11
#define BLOCK_ELINK_MASK           0x07FF
#define BLOCK_ELINK_SHIFT          0
#define BLOCK_GBT_MASK             0x07C0
#define BLOCK_GBT_SHIFT            6
#define BLOCK_EGROUP_MASK          0x0038
#define BLOCK_EGROUP_SHIFT         3
#define BLOCK_EPATH_MASK           0x0007
#define BLOCK_EPATH_SHIFT          0
#define BLOCK_SEQNR_MAX           (BLOCK_SEQNR_MASK>>BLOCK_SEQNR_SHIFT)
#define BLOCK_ELINK_MAX            BLOCK_ELINK_MASK

#define BLOCK_TRAILER_BYTES        2
#define BLOCK_TRAILER_LENGTH_MASK  0x03FF
#define BLOCK_TRAILER_ERROR_MASK   0x0800
#define BLOCK_TRAILER_TRUNC_MASK   0x1000
#define BLOCK_TRAILER_TYPE_MASK    0xE000
#define BLOCK_TRAILER_ERROR_SHIFT  11
#define BLOCK_TRAILER_TRUNC_SHIFT  12
#define BLOCK_TRAILER_TYPE_SHIFT   13


// Chunk descriptor
typedef struct chunk_desc
{
  unsigned int index, length, type;
  bool trunc_or_err, invalid_sz;
} chunk_desc_t;


enum {
  CHUNKTYPE_NULLFILL = 0,
  CHUNKTYPE_FIRST,     // 1->0x20
  CHUNKTYPE_LAST,      // 2->0x40
  CHUNKTYPE_BOTH,      // 3->0x60
  CHUNKTYPE_MIDDLE,    // 4->0x80
  CHUNKTYPE_TIMEOUT,   // 5->0xA0
  CHUNKTYPE_UNDEFINED6,// 6->0xC0
  CHUNKTYPE_OUTOFBAND  // 7->0xE0
};

// Read a file in chunks of 8 MBytes
const int BUFSIZE = 8*1024*1024;

// Matrix template
template <class T, size_t ROW, size_t COL>
using Matrix = std::array<std::array<T, COL>, ROW>;

// Processing matrix. Basically the analysis part with writing to file
template <class T>
std::pair<unsigned, unsigned> processMatrix(T & t, size_t rows, size_t columns)
{
  std::pair<unsigned, unsigned> minmax(0, 0);
  std::ofstream matrixF ("fei4-matrix.txt");
  if (matrixF.is_open())
  {
    float xbar = 0;
    for(size_t i = 0;i < rows; ++i)
    {
      for(size_t j = 0;j < columns; ++j) {
        unsigned val = static_cast<unsigned>(t[i][j]);
        matrixF << val << " ";
        xbar += val;
        if ( val > minmax.second ) {
          minmax.second = val;
        } else if ( val < minmax.first ) {
          minmax.first = val;
        } 
      }
    matrixF << "\n";
    }
    xbar = xbar / 26880; // 80*336
    std::cout << " The mean (avg) is : " << xbar << std::endl;
    
    // find the std dev:
    float numerator = 0;
    float denominator = 26880;
    for (size_t i = 0; i<rows; ++i) {
      for(size_t j = 0;j < columns; ++j) {
        unsigned val = static_cast<float>(t[i][j]);
        numerator = numerator + pow((val - xbar ), 2);
      }
    }


    float standard_deviaton = sqrt( numerator / denominator);
    std::cout << "The standard deviator for the given data is: "
         << standard_deviaton << std::endl;

  } 
  else std::cout << "Unable to open file\n";

  return minmax;
}

// ----------------------------------------------------------------------------



/* This piece of code is highly based on the FlxDataChecker */
unsigned queueInChunks( uint8_t *block, int block_nr,
                        int elink_filter, bool full_mode,
                        std::vector<uint8_t*>& chunks )
{
  //u_short *block_s = (u_short *) block;
  //int elinkseqnr = (int) block_s[0];
  //int gbt   = (elinkseqnr & BLOCK_GBT_MASK)    >> BLOCK_GBT_SHIFT;
  //int grp   = (elinkseqnr & BLOCK_EGROUP_MASK) >> BLOCK_EGROUP_SHIFT;
  //int epath = (elinkseqnr & BLOCK_EPATH_MASK)  >> BLOCK_EPATH_SHIFT;
  //int seqnr = (elinkseqnr & BLOCK_SEQNR_MASK)  >> BLOCK_SEQNR_SHIFT;

  // Apply Elink filter on request
  //  if( elink_filter >= 0 && (elinkseqnr & BLOCK_ELINK_MASK) != elink_filter )
  //    return 0;
  // Go through the chunks in the block, from *end-of-block* to *begin*,
  // collecting info about each for subsequent display from *start-of-block*
  
  uint32_t trailer;
  std::vector<chunk_desc_t> chunk_descs;
  chunk_desc_t chnk;
  int index = BLOCK_BYTES - BLOCK_TRAILER_BYTES;
  while( index > BLOCK_HEADER_BYTES-1 )
    {
      trailer     = block[index] | (block[index+1] << 8);
      chnk.length = trailer & BLOCK_TRAILER_LENGTH_MASK;
      if( full_mode ) chnk.length *= 4;
//      chnk.type   = ((trailer & BLOCK_TRAILER_TYPE_MASK) >>
//                     BLOCK_TRAILER_TYPE_SHIFT);
//      chnk.trunc_or_err = ((trailer & BLOCK_TRAILER_TRUNC_MASK) != 0 ||
//                           (trailer & BLOCK_TRAILER_ERROR_MASK) != 0);

      // Out-Of-Band or Null chunk trailer implies: no payload data
      if( chnk.type == CHUNKTYPE_OUTOFBAND || chnk.type == CHUNKTYPE_NULLFILL )
        {
          chnk.length = 0;
          chnk.trunc_or_err = false;
        }
      // The start of this chunk; account for a possible padding byte
      index -= (chnk.length + (chnk.length & 1));
      chnk.index = index;

      // Is resulting index valid ?
      if( index > BLOCK_HEADER_BYTES-1 )
        {
          chnk.invalid_sz = false;
        }
      else
        {
          // Length can't be correct
          chnk.invalid_sz = true;
          // Adjust for display...
          chnk.index = BLOCK_HEADER_BYTES;
          chnk.length -= (BLOCK_HEADER_BYTES - (index + (chnk.length & 1)));
        }

      //chunks.emplace_back( std::ref(chnk) );
      chunk_descs.push_back( chnk );

      // Move to the preceeding trailer
      index -= BLOCK_TRAILER_BYTES;
   }

  // Collect up chunk pointers.
  int c;
  unsigned int i;
  uint8_t *byt; 
  for( c=chunk_descs.size()-1; c>=0; --c ) {
    byt = &block[chunk_descs[c].index]; 
    for( i=0; i<chunk_descs[c].length; ++i, ++byt ) {
      chunks.push_back( byt );
     }
  }
  return 1;
}


int main( int argc, char *argv[] )
{
  unsigned  blocks_processed  = 0;
  int    elink_filter      = -1;
  bool   full_mode         = false;
  std::string filename;
 
  FILE  *fp;
  unsigned char *buffer = new unsigned char[BUFSIZE];
 
  // Get the file name
  if( optind < argc ) {
    filename = std::string( argv[optind] );
  } else {
    std::cout << "### Filename missing.." << std::endl;
    exit( 0 );
  }

  // Open the file
  if( (fp = fopen( filename.c_str(), "r" )) == NULL ) {
    std::cout << "### Failed to open file " << filename << std::endl;
    exit( 0 );
  }

  // Apply elink filter?
  if( elink_filter >= 0 ) {
    std::cout << "*ELINK FILTER* = "
              << std::hex << std::uppercase << elink_filter << std::dec << std::endl;
  }

  // Get chunks from blocks. 
  std::vector<uint8_t*> chunks;
  unsigned int  blocknr_offs   = 0;
  unsigned long total_bytes    = 0;
  size_t        nbytes;
  while( (nbytes = fread( buffer, 1, BUFSIZE, fp )) > 0 ) {
    for( unsigned int blocknr=0; blocknr<nbytes/BLOCK_BYTES; ++blocknr ) {
       unsigned char *block = &buffer[blocknr*BLOCK_BYTES];
       blocks_processed +=
         queueInChunks( block, blocknr + blocknr_offs,
	                elink_filter, full_mode, chunks );     
    }
    total_bytes  += nbytes;
    blocknr_offs += nbytes/BLOCK_BYTES;
  }
  std::cout << "Blocks processed: " << blocks_processed << std::endl;

  // Data format magic.
  std::vector<uint8_t>  lv1id, bcid, addr, value, 
                             ser_code, ser_num, column, row,
                             tot, cycles_array;  
  Matrix<uint8_t, 80, 336> arrayb, arrayt;
  unsigned lv1id_new = 1000;
  unsigned cycles = 0;
  
  for (unsigned cidx = 0; cidx < chunks.size()/3; ++cidx) {
    uint8_t d1 = (*chunks[cidx*3+1])*256 + (*chunks[cidx*3+2]);   
    if ( *chunks[cidx*3] == 0xe9 ) {
        unsigned ecid_num = (*chunks[cidx*3+1]) / 4;
	unsigned bcid_num = (*chunks[cidx*3+2]) + (static_cast<unsigned>(*chunks[cidx*3+1]) % 4) * 256; 
 
        lv1id.push_back( ecid_num );
        bcid.push_back( bcid_num );
        
        if ( lv1id_new != ecid_num ) { cycles = 0; }
        else { cycles++; }

        lv1id_new = ecid_num;
        lv1id.push_back( (*chunks[cidx*3+1])>>2 );
 
    } else if ( *chunks[cidx*3] == 0xea ) { addr.push_back( d1 );
    } else if ( *chunks[cidx*3] == 0xec ) { value.push_back( d1 ); 
    } else if ( *chunks[cidx*3] == 0xef ) {
      ser_code.push_back( (*chunks[cidx*3+1]) / 4 );
      ser_num.push_back( (*chunks[cidx*3+2]) + (static_cast<unsigned>(*chunks[cidx*3+1]) % 4) * 256 );
    } else if ( *chunks[cidx*3] == 0xab ) {
      std::cout << "Empty Record is found! \n  -> Empty record should not exist in 8b10b mode!\n"; 
      break;
    } else {

      unsigned col_dat = (*chunks[cidx*3])/2;
      unsigned row_dat = ((*chunks[cidx*3])%2)*256 + (*chunks[cidx*3+1]);
      unsigned tot1 = (*chunks[cidx*3+2])/16;
      unsigned tot2 = (*chunks[cidx*3+2])%16;

      uint8_t tot1_real, tot2_real;
      if ( tot2 == 15 ) {
        if ( tot1 == 14 ){ 
          tot1_real = 0; 
        } else if ( tot1 == 15 ) {
          std::cout << "TOT data error!\n";                    
          exit(1);
        } else {
          tot1_real = tot1 + 1;
        }

        // APPENDS
        column.push_back( col_dat );
        row.push_back( row_dat );
        tot.push_back( tot1_real );
        cycles_array.push_back( cycles );
        arrayb[col_dat-1][row_dat-1] += 1;
        arrayt[col_dat-1][row_dat-1] += tot1_real;
        // if (...) print to file is missing from python ana.

      } else { // tot2!=15
      
        if ( tot1 == 14 ){
          tot1_real = 0;
        } else if ( tot1 == 15 ) {
          std::cout << "TOT data error!\n";
          exit(1);
        } else {
          tot1_real = tot1 + 1;
        }

        if ( tot2 == 14 ) { tot2_real = 0; }
        else { tot2_real = tot2 + 1; }

        column.push_back( col_dat );
        column.push_back( col_dat );
        row.push_back( row_dat );
        row.push_back( row_dat + 1 );
        tot.push_back( tot1_real );
        tot.push_back( tot2_real );
        cycles_array.push_back( cycles );
        cycles_array.push_back( cycles );
        arrayb[col_dat-1][row_dat-1] += 1;
        arrayb[col_dat-1][row_dat] += 1;
        arrayt[col_dat-1][row_dat-1] += tot1_real;
        arrayt[col_dat-1][row_dat] += tot2_real;

      }
    }
  }

  // Process matrix and print min/max triggered. 
  std::pair<unsigned, unsigned> minmax = processMatrix( arrayb, 80, 336 );
  std::cout << "Min triggered: " << minmax.first << " | Max triggered: " << minmax.second;

}

// ----------------------------------------------------------------------------


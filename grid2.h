// grid.h.
// Organism class using a simple grid of tiles

#ifndef ORGANISM_H
#define ORGANISM_H

#include <cstddef>
#include <cstdint>
#include <vector>
#include <map>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include <cstring>
#include <iostream>

// Uncomment next line to use Howard Hinnant Short Allocator
// #define USE_HOWARDH_ALLOC_L 1

#ifdef USE_HOWARDH_ALLOC_L
#include "short_alloc.h"
#endif

// Make GCC __builtin_popcountll available with MSVC
// (used to calculate the Hamming weight efficiently)
#ifdef _MSC_VER
#include <intrin.h>
#define __builtin_popcountll __popcnt64
#define __builtin_popcount   __popcnt
#endif

static_assert(sizeof(size_t) == sizeof(int64_t), "size_t too small, need a 64-bit compile");

// CELL
typedef int32_t  cell_coord_type;
typedef uint64_t cell_whole_type;
typedef int64_t  cell_signed_whole_type;

#if 0
#define XX(w)    ((cell_coord_type)(w))
#define YY(w)    ((cell_coord_type)(((w) >> 32) & 0xFFFFFFFF))
#define WW(x, y) (((cell_signed_whole_type)(y) << 32) | ((cell_signed_whole_type)(x) & 0xFFFFFFFF))
#else
#define XX(w)    ((cell_coord_type)(((w) >> 32) & 0xFFFFFFFF))
#define YY(w)    ((cell_coord_type)(w))
#define WW(x, y) (((cell_signed_whole_type)(x) << 32) | ((cell_signed_whole_type)(y) & 0xFFFFFFFF))
#endif

struct Cell {
   cell_coord_type x;
   cell_coord_type y;
   Cell(cell_coord_type x_, cell_coord_type y_) : x(x_), y(y_) {}
   Cell(cell_whole_type w) : x(XX(w)), y(YY(w)) {}
};
inline bool operator<(const Cell& lhs, const Cell& rhs) {
   if (lhs.x < rhs.x) return true;
   if (lhs.x > rhs.x) return false;
   return lhs.y < rhs.y;
   // ... or one-liner: return lhs.x < rhs.x || (lhs.x == rhs.x && lhs.y < rhs.y);
}
inline bool operator>( const Cell& lhs, const Cell& rhs) { return rhs < lhs; }
inline bool operator<=(const Cell& lhs, const Cell& rhs) { return !(lhs > rhs); }
inline bool operator>=(const Cell& lhs, const Cell& rhs) { return !(lhs < rhs); }
inline bool operator==(const Cell& lhs, const Cell& rhs) { return lhs.x == rhs.x && lhs.y == rhs.y; }
inline bool operator!=(const Cell& lhs, const Cell& rhs) { return !(lhs == rhs); }
using cell_list_type = std::vector<Cell>;

// -----------------------------------------------------------
// SQUARE TILE
// Set TILE_SIZE_FULL to either 32 or 64
#define TILE_SIZE_FULL    64

#if TILE_SIZE_FULL == 64
typedef uint64_t uintsq_t;
#define BM_POPCOUNT         __builtin_popcountll
#define BM_MIDDLE           0x3ffffffffffffffcull
#define BM_LEFT             0xfffffffffffffffcull
#define BM_RIGHT            0x3fffffffffffffffull
#define BM_OUTER            0xc000000000000003ull
#define BM_LEFTMIDDLE       0x3000000000000000ull
#define BM_RIGHTMIDDLE      0x000000000000000cull
#elif TILE_SIZE_FULL == 32
typedef uint32_t uintsq_t;
#define BM_POPCOUNT         __builtin_popcount
#define BM_MIDDLE           0x3ffffffc
#define BM_LEFT             0xfffffffc
#define BM_RIGHT            0x3fffffff
#define BM_OUTER            0xc0000003
#define BM_LEFTMIDDLE       0x30000000
#define BM_RIGHTMIDDLE      0x0000000c
#else
#error "invalid TILE_SIZE_FULL"
#endif

#define BORDER_WIDTH        2
#define BORDER_WIDTH_P1     (BORDER_WIDTH + 1)
#define TILE_SIZE_FULL_M1   (TILE_SIZE_FULL - 1)
#define TILE_SIZE_MBD       (TILE_SIZE_FULL - BORDER_WIDTH)
#define TILE_SIZE_MBD_M1    (TILE_SIZE_MBD - 1)
#define TILE_SIZE_CORE      (TILE_SIZE_FULL - 2 * BORDER_WIDTH)
#define TILE_SIZE_CORE_P1   (TILE_SIZE_CORE + 1)

// Neighbours are numbered clockwise starting with the one directly above
#define NUM_NEIGH            8
#define NEIGH_TOP            0
#define NEIGH_TOP_RIGHT      1
#define NEIGH_RIGHT          2
#define NEIGH_BOTTOM_RIGHT   3
#define NEIGH_BOTTOM         4
#define NEIGH_BOTTOM_LEFT    5
#define NEIGH_LEFT           6
#define NEIGH_TOP_LEFT       7
#define NEIGH_ANY            0xff
#define NEIGH_BIT(i)         (1 << (i))
#define IS_NEIGH(n, i)       ((n) & NEIGH_BIT(i))
// Note: i ^ 4 trick sets NEIGH flag to the opposite one:
//   0 > 4  (TOP > BOTTOM)
//   1 > 5  (TOP RIGHT > BOTTOM LEFT)
//   2 > 6  (RIGHT > LEFT)
//   3 > 7  (BOTTOM RIGHT > TOP LEFT)
//   4 > 0  (BOTTOM > TOP)
//   5 > 1  (BOTTOM LEFT > TOP RIGHT)
//   6 > 2  (LEFT > RIGHT)
//   7 > 3  (TOP LEFT > BOTTOM RIGHT)

// Note: mapping of x (cell) to tx (tile) is:
//        x       tx
//   ----------   --
//   ...
//   -121..-180   -3
//    -61..-120   -2
//     -1..-60    -1
//      0.. 59     0
//     60..119     1
//    120..179     2
//   ...
// Ditto for y (cell) to ty (tile).

// Converse of get_tile_coords()
// Given tx/ty and ix/iy return x/y the cell coordinate
#define GET_CELL_COORD(tx, ix) (TILE_SIZE_CORE * (tx) + (ix) - BORDER_WIDTH)

// Input cell (x, y). Return (tx, ty, ix, iy)
// (tx, ty) : Tile coords
// (ix, iy) : x, y coords inside tile
static void get_tile_coords(
   cell_coord_type x, cell_coord_type y,
   cell_coord_type& tx, cell_coord_type& ty, cell_coord_type& ix, cell_coord_type& iy
) {
   cell_coord_type ox = x % TILE_SIZE_CORE;
   if (ox < 0) ox += TILE_SIZE_CORE;
   tx = (x - ox) / TILE_SIZE_CORE;
   ix = ox + BORDER_WIDTH;

   cell_coord_type oy = y % TILE_SIZE_CORE;
   if (oy < 0) oy += TILE_SIZE_CORE;
   ty = (y - oy) / TILE_SIZE_CORE;
   iy = oy + BORDER_WIDTH;
}

// A simple square tile bitmap (64 x 64 or 32 x 32)
// x and y must be in 0..TILE_SIZE_FULL_M1 range
// Cell value in row[] bitmap is 0 (dead) or 1 (alive)
struct SqTile {
   uintsq_t row[TILE_SIZE_FULL];
   SqTile() { memset(row, 0, sizeof(row)); }
   // bool is_empty(int& top, int& bottom) const {
   //    top = 0;
   //    bottom = TILE_SIZE_FULL_M1;
   //    while (top < TILE_SIZE_FULL_M1 && row[top] == 0) ++top;
   //    while (bottom > 0 && row[bottom] == 0) --bottom;
   //    return top > bottom;
   // }
   int getcellval(int x, int y) const {
      uintsq_t mask = (uintsq_t)1 << (TILE_SIZE_FULL_M1 - x);
      return row[y] & mask ? 1 : 0;
   }
   void setcellval(int x, int y, int v) {
      uintsq_t mask = (uintsq_t)1 << (TILE_SIZE_FULL_M1 - x);
      if (v) { row[y] |= mask; } else { row[y] &= ~mask; }
   }
   // Used to initialize starting state
   void insertcells(const cell_list_type& cells) {
      for (const auto& c : cells) setcellval(c.x, c.y, 1);
   }
   // Used for verification and testing
   uintsq_t count() const {
      uintsq_t cnt = 0;
      for (int y = 0; y < TILE_SIZE_FULL; ++y) {
         if (!row[y]) continue;
         cnt += BM_POPCOUNT(row[y]);
      }
      return cnt;
   }
   // Used for verification and testing
   cell_list_type getlivecells() const
   {
      cell_list_type cells;
      for (int y = 0; y < TILE_SIZE_FULL; ++y) {
         if (!row[y]) continue;
         for (int x = 0; x < TILE_SIZE_FULL; ++x) {
            if (getcellval(x, y)) { cells.push_back(Cell(x, y)); }
         }
      }
      std::sort(cells.begin(), cells.end());
      return cells;
   }

   // Advance the interior of square tile by one tick.
   // Return 1 if square tile changed, else 0.
   // neigh is an out parameter containing update flags
   // of the eight neighbours of this square tile.
   int tiletick(int top, int& neigh) {
      neigh = 0;
      uintsq_t bigdiff = 0;
      uintsq_t carry[TILE_SIZE_FULL];
      uintsq_t parity[TILE_SIZE_FULL];
      uintsq_t diff[TILE_SIZE_FULL];
      memset(carry,  0, sizeof(carry));
      memset(parity, 0, sizeof(parity));
      memset(diff,   0, sizeof(diff));
      uintsq_t aa, bb, p, q, r, s, bit0, bit1, bit2;
#if 1
      // int top = 0;
      // while (top < TILE_SIZE_FULL_M1 && row[top] == 0) ++top;
      int bottom = TILE_SIZE_FULL_M1;
      while (bottom > 0 && row[bottom] == 0) --bottom;
      // if (top > bottom) return 0;   // cannot happen with get_top check
#endif

      for (int i = top; i <= bottom; ++i) {
         aa = row[i] >> 1;
         bb = row[i] << 1;
         q = aa ^ bb;
         parity[i] = q ^ row[i];
         carry[i] = (q & row[i]) | (aa & bb);
      }
      --top; ++bottom;
      if (top < 1) top = 1;
      if (bottom > TILE_SIZE_MBD) bottom = TILE_SIZE_MBD;
      for (int i = top; i <= bottom; ++i) {
         aa = parity[i-1];
         bb = parity[i+1];
         q = aa ^ bb;
         bit0 = q ^ parity[i];
         r = (q & parity[i]) | (aa & bb);

         aa = carry[i-1];
         bb = carry[i+1];
         q = aa ^ bb;
         p = q ^ carry[i];
         s = (q & carry[i]) | (aa & bb);

         bit1 = p ^ r;
         bit2 = s ^ (p & r);
         p = (bit0 & bit1 & ~bit2) | (bit2 & ~bit1 & ~bit0 & row[i]);
         diff[i] = (row[i] ^ p) & BM_MIDDLE;
         bigdiff |= diff[i];
         row[i] = (p & BM_MIDDLE) | (row[i] & ~BM_MIDDLE);
      }
      aa = diff[BORDER_WIDTH]   | diff[BORDER_WIDTH_P1];
      bb = diff[TILE_SIZE_CORE] | diff[TILE_SIZE_CORE_P1];
      if (bigdiff) {
         if (bigdiff & BM_LEFTMIDDLE)  neigh |= NEIGH_BIT(NEIGH_LEFT);
         if (bigdiff & BM_RIGHTMIDDLE) neigh |= NEIGH_BIT(NEIGH_RIGHT);
      }
      if (aa) {
         neigh |= NEIGH_BIT(NEIGH_TOP);
         if (aa & BM_LEFTMIDDLE)  neigh |= NEIGH_BIT(NEIGH_TOP_LEFT);
         if (aa & BM_RIGHTMIDDLE) neigh |= NEIGH_BIT(NEIGH_TOP_RIGHT);
      }
      if (bb) {
         neigh |= NEIGH_BIT(NEIGH_BOTTOM);
         if (bb & BM_LEFTMIDDLE)  neigh |= NEIGH_BIT(NEIGH_BOTTOM_LEFT);
         if (bb & BM_RIGHTMIDDLE) neigh |= NEIGH_BIT(NEIGH_BOTTOM_RIGHT);
      }
      return (bigdiff != 0) ? 1 : 0;
   }

   // Advance the interior of square tile by two ticks.
   // Return 1 if square tile changed, else 0.
   // neigh is an out parameter containing update flags
   // of the eight neighbours of this square tile.
   int tiletwoticks(int top, int& neigh) {
      neigh = 0;
      uintsq_t bigdiff = 0;
      uintsq_t carry[TILE_SIZE_FULL];
      uintsq_t parity[TILE_SIZE_FULL];
      uintsq_t diff[TILE_SIZE_FULL];
      uintsq_t ee[TILE_SIZE_FULL];
      memset(carry,  0, sizeof(carry));
      memset(parity, 0, sizeof(parity));
      memset(diff,   0, sizeof(diff));
      memset(ee,     0, sizeof(ee));
      uintsq_t aa, bb, p, q, r, s, bit0, bit1, bit2;
#if 1
      // int top = 0;
      // while (top < TILE_SIZE_FULL_M1 && row[top] == 0) ++top;
      int bottom = TILE_SIZE_FULL_M1;
      while (bottom > 0 && row[bottom] == 0) --bottom;
      // if (top > bottom) return 0;   // cannot happen with get_top check
#endif

      // even --> odd
      for (int i = top; i <= bottom; ++i) {
         aa = row[i] >> 1;
         bb = row[i] << 1;
         q = aa ^ bb;
         parity[i] = q ^ row[i];
         carry[i] = (q & row[i]) | (aa & bb);
      }
      --top; ++bottom;
      if (top < 1) top = 1;
      if (bottom > TILE_SIZE_MBD) bottom = TILE_SIZE_MBD;

      for (int i = top; i <= bottom; ++i) {
         aa = parity[i-1];
         bb = parity[i+1];
         q = aa ^ bb;
         bit0 = q ^ parity[i];
         r = (q & parity[i]) | (aa & bb);

         aa = carry[i-1];
         bb = carry[i+1];
         q = aa ^ bb;
         p = q ^ carry[i];
         s = (q & carry[i]) | (aa & bb);

         bit1 = p ^ r;
         bit2 = s ^ (p & r);
         ee[i] = (bit0 & bit1 & ~bit2) | (bit2 & ~bit1 & ~bit0 & row[i]);
      }

      // odd --> even
      for (int i = top; i <= bottom; ++i) {
         aa = ee[i] >> 1;
         bb = ee[i] << 1;
         q = aa ^ bb;
         parity[i] = q ^ ee[i];
         carry[i] = (q & ee[i]) | (aa & bb);
      }
      --top; ++bottom;
      if (top < BORDER_WIDTH) top = BORDER_WIDTH;
      if (bottom > TILE_SIZE_MBD_M1) bottom = TILE_SIZE_MBD_M1;

      for (int i = top; i <= bottom; ++i) {
         aa = parity[i-1];
         bb = parity[i+1];
         q = aa ^ bb;
         bit0 = q ^ parity[i];
         r = (q & parity[i]) | (aa & bb);

         aa = carry[i-1];
         bb = carry[i+1];
         q = aa ^ bb;
         p = q ^ carry[i];
         s = (q & carry[i]) | (aa & bb);

         bit1 = p ^ r;
         bit2 = s ^ (p & r);
         p = (bit0 & bit1 & ~bit2) | (bit2 & ~bit1 & ~bit0 & ee[i]);
         diff[i] = (row[i] ^ p) & BM_MIDDLE;
         bigdiff |= diff[i];
         row[i] = (p & BM_MIDDLE) | (row[i] & ~BM_MIDDLE);
      }

      aa = diff[BORDER_WIDTH]   | diff[BORDER_WIDTH_P1];
      bb = diff[TILE_SIZE_CORE] | diff[TILE_SIZE_CORE_P1];
      if (bigdiff) {
         if (bigdiff & BM_LEFTMIDDLE)  neigh |= NEIGH_BIT(NEIGH_LEFT);
         if (bigdiff & BM_RIGHTMIDDLE) neigh |= NEIGH_BIT(NEIGH_RIGHT);
      }
      if (aa) {
         neigh |= NEIGH_BIT(NEIGH_TOP);
         if (aa & BM_LEFTMIDDLE)  neigh |= NEIGH_BIT(NEIGH_TOP_LEFT);
         if (aa & BM_RIGHTMIDDLE) neigh |= NEIGH_BIT(NEIGH_TOP_RIGHT);
      }
      if (bb) {
         neigh |= NEIGH_BIT(NEIGH_BOTTOM);
         if (bb & BM_LEFTMIDDLE)  neigh |= NEIGH_BIT(NEIGH_BOTTOM_LEFT);
         if (bb & BM_RIGHTMIDDLE) neigh |= NEIGH_BIT(NEIGH_BOTTOM_RIGHT);
      }
      return (bigdiff != 0) ? 1 : 0;
   }
};

// A more complex square tile.
struct SquareTile {
   SquareTile() = default;
   SqTile sq;
   cell_coord_type tx;
   cell_coord_type ty;
   int updateflags;
   struct SquareTile* neighbours[NUM_NEIGH];
};

#ifdef USE_HOWARDH_ALLOC_L
// Create a vector<T> template with a small buffer of 512 bytes.
//   Note for vector it is possible to reduce the alignment requirements
//   down to alignof(T) because vector doesn't allocate anything but T's.
//   And if we're wrong about that guess, it is a compile-time error, not
//   a run time error.
template <class T, std::size_t BufSize = 512>
using SmallVector = std::vector<T, short_alloc<T, BufSize, alignof(T)>>;
#endif

struct TileHasher { size_t operator()(cell_whole_type w) const { return w; }  };
using tile_map_type = std::unordered_map<cell_whole_type, SquareTile, TileHasher>;
using tile_sorted_map_type = std::map<cell_whole_type, SquareTile>;
using tile_list_type = std::vector<SquareTile*>;
#ifdef USE_HOWARDH_ALLOC_L
using tile_list_arena_type = SmallVector<SquareTile*>;
#endif

static void dump_one_tile(const SquareTile& sqt) {
   size_t popc = 0;
   for (int iy = BORDER_WIDTH; iy <= TILE_SIZE_MBD_M1; ++iy) {
      if (!sqt.sq.row[iy]) continue;
      popc += BM_POPCOUNT(sqt.sq.row[iy] & BM_MIDDLE);
   }
   std::cerr << "tx=" << sqt.tx << " ty=" << sqt.ty << " popc=" << popc << "\n";
   std::cerr << "  updateflags=";
   for (int i = 31; i >= 0; --i) std::cerr << ((sqt.updateflags >> i) & 1);
   std::cerr << "\n";
   std::cerr << "  live neighbours:";
   if (sqt.neighbours) {
      for (int n = 0; n < NUM_NEIGH; ++n) {
         if (!sqt.neighbours[n]) continue;
         std::cerr << " " << n;
      }
   }
   std::cerr << "\n";
}

// -----------------------------------------------------------
// ORGANISM
class Organism {
public:
   // Organism(size_t ntiles = 2048) {
   //    tiles_m.reserve(ntiles);
   //    modified_m.reserve(64);
   // }

   // Used to initialize the starting state of the organism
   void insert_cells(const cell_list_type& cells) {
      for (const auto& c : cells) setcell(c.x, c.y, 1);
   }

   // Used for verification and testing the state of the organism
   void dump_tiles() const {
      tile_sorted_map_type ordered(tiles_m.cbegin(), tiles_m.cend());
      size_t ntiles = tiles_m.size();
      std::cerr << "=== Dump Tiles, n=" << ntiles << " ==========\n";
      size_t ii = 0;
      // for (const auto& tt : tiles_m) {
      for (const auto& tt : ordered) {
         ++ii;
         // cell_whole_type k = tt.first;
         // Note: kx,ky match Tx,Ty
         // cell_coord_type kx = XX(k);
         // cell_coord_type ky = YY(k);
         std::cerr << ii << ":-----------------------------\n";
         const SquareTile& sqt = tt.second;
         dump_one_tile(sqt);
      }
      size_t nmodified = modified_m.size();
      std::cerr << "=== Dump Modif, n=" << nmodified << " ==========\n";
      ii = 0;
      for (auto sqt : modified_m) {
         ++ii;
         std::cerr << ii << ":-----------------------------\n";
         dump_one_tile(*sqt);
      }
   }

   // Used for verification and testing the state of the organism
   size_t count() const {
      size_t cnt = 0;
      for (const auto& tt : tiles_m) {
         const SquareTile& sqt = tt.second;
         for (int iy = BORDER_WIDTH; iy <= TILE_SIZE_MBD_M1; ++iy) {
            if (!sqt.sq.row[iy]) continue;
            cnt += BM_POPCOUNT(sqt.sq.row[iy] & BM_MIDDLE);
         }
      }
      return cnt;
   }

   // Used for verification and testing the state of the organism
   cell_list_type get_live_cells() const
   {
      cell_list_type cells;
      for (const auto& tt : tiles_m) {
         const SquareTile& sqt = tt.second;
         for (int iy = BORDER_WIDTH; iy <= TILE_SIZE_MBD_M1; ++iy) {
            if (!sqt.sq.row[iy]) continue;
            for (int ix = BORDER_WIDTH; ix <= TILE_SIZE_MBD_M1; ++ix) {
               if (sqt.sq.getcellval(ix, iy)) {
                  cells.push_back(Cell( GET_CELL_COORD(sqt.tx, ix),
                                        GET_CELL_COORD(sqt.ty, iy) ));
               }
            }
         }
      }
      std::sort(cells.begin(), cells.end());
      return cells;
   }

   SquareTile* get_neighbour(SquareTile* sqt, int i) {
      if (!sqt->neighbours[i]) {
         cell_coord_type x = sqt->tx;
         cell_coord_type y = sqt->ty;
         if (i >= NEIGH_TOP_RIGHT    && i <= NEIGH_BOTTOM_RIGHT) ++x;
         if (i >= NEIGH_BOTTOM_RIGHT && i <= NEIGH_BOTTOM_LEFT)  ++y;
         if (i >= NEIGH_BOTTOM_LEFT  && i <= NEIGH_TOP_LEFT)     --x;
         if (i == NEIGH_TOP_LEFT     || i <= NEIGH_TOP_RIGHT)    --y;
         // Note: next line creates a new entry in tiles_m[]
         // if the key does not already exist
         sqt->neighbours[i] = &tiles_m[WW(x, y)];
         sqt->neighbours[i]->tx = x;
         sqt->neighbours[i]->ty = y;
      }
      return sqt->neighbours[i];
   }

   // Alert the neighbour that its neighbour (the original tile) has changed
   void update_neighbour(SquareTile* sqt, int i) {
      if (get_neighbour(sqt, i)->updateflags == 0)
         modified_m.push_back(get_neighbour(sqt, i));
      get_neighbour(sqt, i)->updateflags |= NEIGH_BIT(i ^ 4);
   }

   // XXX: update_tile and update_tile_two are identical except
   // for calling tiletick() versus tiletwoticks()
   void update_tile(SquareTile* sqt) {
      int neigh, update_flag;
      int top = 0;
      while (top < TILE_SIZE_FULL && sqt->sq.row[top] == 0) ++top;
      if (top == TILE_SIZE_FULL) return;
      update_flag = sqt->sq.tiletick(top, neigh);
      if (update_flag) {
         if (sqt->updateflags == 0) modified_m.push_back(sqt);
         sqt->updateflags |= NEIGH_BIT(NUM_NEIGH);
      }
      if (!neigh) return;
      for (int i = 0; i < NUM_NEIGH; ++i) {
         if (neigh & NEIGH_BIT(i)) update_neighbour(sqt, i);
      }
   }

   void update_tile_two(SquareTile* sqt) {
      int neigh, update_flag;
      int top = 0;
      while (top < TILE_SIZE_FULL && sqt->sq.row[top] == 0) ++top;
      if (top == TILE_SIZE_FULL) return;
      update_flag = sqt->sq.tiletwoticks(top, neigh);
      if (update_flag) {
         if (sqt->updateflags == 0) modified_m.push_back(sqt);
         sqt->updateflags |= NEIGH_BIT(NUM_NEIGH);
      }
      if (!neigh) return;
      for (int i = 0; i < NUM_NEIGH; ++i) {
         if (neigh & NEIGH_BIT(i)) update_neighbour(sqt, i);
      }
   }

   // Update the relevant portions of the boundary (a 64-by-64 square
   // with the central 60-by-60 square removed) by copying data from
   // the interiors (the 60-by-60 central squares) of the neighbours.
   // Only perform this copying when necessary.
   // Note: alternatively: 32-by-32 with central 28-by-28.
   void update_boundary_tile(SquareTile* sqt) {
      if (sqt->updateflags & NEIGH_BIT(NEIGH_TOP)) {
         SquareTile* n = get_neighbour(sqt, NEIGH_TOP);
         sqt->sq.row[0] = (n->sq.row[TILE_SIZE_CORE] & BM_MIDDLE) | (sqt->sq.row[0] & BM_OUTER);
         sqt->sq.row[1] = (n->sq.row[TILE_SIZE_CORE_P1] & BM_MIDDLE) | (sqt->sq.row[1] & BM_OUTER);
      }
      if (sqt->updateflags & NEIGH_BIT(NEIGH_TOP_LEFT)) {
         SquareTile* n = get_neighbour(sqt, NEIGH_TOP_LEFT);
         sqt->sq.row[0] = ((n->sq.row[TILE_SIZE_CORE] & BM_MIDDLE) << TILE_SIZE_CORE) | (sqt->sq.row[0] & BM_RIGHT);
         sqt->sq.row[1] = ((n->sq.row[TILE_SIZE_CORE_P1] & BM_MIDDLE) << TILE_SIZE_CORE) | (sqt->sq.row[1] & BM_RIGHT);
      }
      if (sqt->updateflags & NEIGH_BIT(NEIGH_TOP_RIGHT)) {
         SquareTile* n = get_neighbour(sqt, NEIGH_TOP_RIGHT);
         sqt->sq.row[0] = ((n->sq.row[TILE_SIZE_CORE] & BM_MIDDLE) >> TILE_SIZE_CORE) | (sqt->sq.row[0] & BM_LEFT);
         sqt->sq.row[1] = ((n->sq.row[TILE_SIZE_CORE_P1] & BM_MIDDLE) >> TILE_SIZE_CORE) | (sqt->sq.row[1] & BM_LEFT);
      }
      if (sqt->updateflags & NEIGH_BIT(NEIGH_BOTTOM)) {
         SquareTile* n = get_neighbour(sqt, NEIGH_BOTTOM);
         sqt->sq.row[TILE_SIZE_MBD] = (n->sq.row[BORDER_WIDTH] & BM_MIDDLE) | (sqt->sq.row[TILE_SIZE_MBD] & BM_OUTER);
         sqt->sq.row[TILE_SIZE_FULL_M1] = (n->sq.row[3] & BM_MIDDLE) | (sqt->sq.row[TILE_SIZE_FULL_M1] & BM_OUTER);
      }
      if (sqt->updateflags & NEIGH_BIT(NEIGH_BOTTOM_LEFT)) {
         SquareTile* n = get_neighbour(sqt, NEIGH_BOTTOM_LEFT);
         sqt->sq.row[TILE_SIZE_MBD] = ((n->sq.row[BORDER_WIDTH] & BM_MIDDLE) << TILE_SIZE_CORE) | (sqt->sq.row[TILE_SIZE_MBD] & BM_RIGHT);
         sqt->sq.row[TILE_SIZE_FULL_M1] = ((n->sq.row[3] & BM_MIDDLE) << TILE_SIZE_CORE) | (sqt->sq.row[TILE_SIZE_FULL_M1] & BM_RIGHT);
      }
      if (sqt->updateflags & NEIGH_BIT(NEIGH_BOTTOM_RIGHT)) {
         SquareTile* n = get_neighbour(sqt, NEIGH_BOTTOM_RIGHT);
         sqt->sq.row[TILE_SIZE_MBD] = ((n->sq.row[BORDER_WIDTH] & BM_MIDDLE) >> TILE_SIZE_CORE) | (sqt->sq.row[TILE_SIZE_MBD] & BM_LEFT);
         sqt->sq.row[TILE_SIZE_FULL_M1] = ((n->sq.row[3] & BM_MIDDLE) >> TILE_SIZE_CORE) | (sqt->sq.row[TILE_SIZE_FULL_M1] & BM_LEFT);
      }
      if (sqt->updateflags & NEIGH_BIT(NEIGH_LEFT)) {
         SquareTile* n = get_neighbour(sqt, NEIGH_LEFT);
         for (int i = BORDER_WIDTH; i < TILE_SIZE_MBD; ++i) {
            sqt->sq.row[i] = ((n->sq.row[i] & BM_MIDDLE) << TILE_SIZE_CORE) | (sqt->sq.row[i] & BM_RIGHT);
         }
      }
      if (sqt->updateflags & NEIGH_BIT(NEIGH_RIGHT)) {
         SquareTile* n = get_neighbour(sqt, NEIGH_RIGHT);
         for (int i = BORDER_WIDTH; i < TILE_SIZE_MBD; ++i) {
            sqt->sq.row[i] = ((n->sq.row[i] & BM_MIDDLE) >> TILE_SIZE_CORE) | (sqt->sq.row[i] & BM_LEFT);
         }
      }
   }

   void tick() {
      // Update boundary of all modified tiles
      for (auto sqt : modified_m) {
         if (sqt->updateflags & NEIGH_ANY) update_boundary_tile(sqt);
         sqt->updateflags = 0;
      }

      // Create old_modified (used to create new modified list)
      // Next two are slower:
      //   tile_list_type old_modified; old_modified.swap(modified_m);
      //   tile_list_type old_modified = std::move(modified_m);
#ifdef USE_HOWARDH_ALLOC_L
      tile_list_arena_type::allocator_type::arena_type st_arena;
      const tile_list_arena_type old_modified(modified_m.cbegin(), modified_m.cend(), st_arena);
#else
      const tile_list_type old_modified(modified_m.cbegin(), modified_m.cend());
#endif
      modified_m.clear();

      // Update core of all modified tiles, creating new modified list
      for (auto sqt : old_modified) { update_tile(sqt); }
   }

   void twoticks() {
      // Update boundary of all modified tiles
      for (auto sqt : modified_m) {
         if (sqt->updateflags & NEIGH_ANY) update_boundary_tile(sqt);
         sqt->updateflags = 0;
      }

      // Create old_modified (used to create new modified list)
      // Next two are slower:
      //   tile_list_type old_modified; old_modified.swap(modified_m);
      //   tile_list_type old_modified = std::move(modified_m);
#ifdef USE_HOWARDH_ALLOC_L
      tile_list_arena_type::allocator_type::arena_type st_arena;
      const tile_list_arena_type old_modified(modified_m.cbegin(), modified_m.cend(), st_arena);
#else
      const tile_list_type old_modified(modified_m.cbegin(), modified_m.cend());
#endif
      modified_m.clear();

      // Update core of all modified tiles, creating new modified list
      for (auto sqt : old_modified) { update_tile_two(sqt); }
   }

   void ticks(size_t nticks) {
      if (nticks == 1) { tick(); return; }
      size_t half = nticks >> 1;
      size_t rem  = nticks &  1;
      for (size_t i = 0; i < half; ++i) { twoticks(); }
      if (rem) { tick(); }
   }

   void updatecell(SquareTile* sqt, cell_coord_type x, cell_coord_type y) {
      if (sqt->updateflags == 0) modified_m.push_back(sqt);
      sqt->updateflags |= NEIGH_BIT(NUM_NEIGH);
      if (y <= BORDER_WIDTH_P1) update_neighbour(sqt, NEIGH_TOP);
      if (y >= TILE_SIZE_CORE) update_neighbour(sqt, NEIGH_BOTTOM);
      if (x <= BORDER_WIDTH_P1) {
         update_neighbour(sqt, NEIGH_LEFT);
         if (y <= BORDER_WIDTH_P1) update_neighbour(sqt, NEIGH_TOP_LEFT);
         if (y >= TILE_SIZE_CORE) update_neighbour(sqt, NEIGH_BOTTOM_LEFT);
      }
      if (x >= TILE_SIZE_CORE) {
         update_neighbour(sqt, NEIGH_RIGHT);
         if (y <= BORDER_WIDTH_P1) update_neighbour(sqt, NEIGH_TOP_RIGHT);
         if (y >= TILE_SIZE_CORE) update_neighbour(sqt, NEIGH_BOTTOM_RIGHT);
      }
   }

   void setcell(cell_coord_type x, cell_coord_type y, int state) {
      cell_coord_type tx, ty, ix, iy;
      get_tile_coords(x, y, tx, ty, ix, iy);
      // Note: next line creates a new entry in tiles_m[]
      // if the key does not already exist
      SquareTile* sqt = &tiles_m[WW(tx, ty)];
      sqt->tx = tx;
      sqt->ty = ty;
      sqt->sq.setcellval(ix, iy, state);
      updatecell(sqt, ix, iy);
   }

   int getcellval(cell_coord_type x, cell_coord_type y) const {
      cell_coord_type tx, ty, ix, iy;
      get_tile_coords(x, y, tx, ty, ix, iy);
      auto kk = tiles_m.find(WW(tx, ty));
      if (kk == tiles_m.end()) { return 0; }
      return kk->second.sq.getcellval(ix, iy);
   }

private:
   tile_map_type tiles_m;        // For lidka up to ~2048 entries
   tile_list_type modified_m;    // For lidka up to ~64 entries
   // Could try USE_HOWARDH_ALLOC_L (but seems to be slower here)
   // tile_list_arena_type::allocator_type::arena_type st_arena_m;
   // tile_list_arena_type modified_m {st_arena_m};
};

#endif

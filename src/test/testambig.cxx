/// Test isodual3D_ambig routines

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2012 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "isodual3D_ambig.h"

#include "ijkstring.txx"

using namespace std;

// global variables

// output routines
void out_ambig_table(const ISODUAL3D::AMBIG_TABLE & ambig_table);

// check routines
void check_ambig_table(const ISODUAL3D::AMBIG_TABLE & ambig_table);

// miscellaneous routines
void memory_exhaustion();
void help(), usage_error();
void parse_command_line(int argc, char **argv);


int main(int argc, char **argv)
{

  try {

    std::set_new_handler(memory_exhaustion);

    parse_command_line(argc, argv);

    ISODUAL3D::AMBIG_TABLE ambig_table;
    ambig_table.SetCubeAmbiguityTable();

    out_ambig_table(ambig_table);

    check_ambig_table(ambig_table);
    cout << "Passed check." << endl;
  }
  catch (IJK::ERROR error) {
    if (error.NumMessages() == 0) {
      cerr << "Unknown error." << endl;
    }
    else { error.Print(cerr); }
    cerr << "Exiting." << endl;
    exit(30);
  }
  catch (...) {
    cerr << "Unknown error." << endl;
    exit(50);
  };

}

// **************************************************
// I/O routines
// **************************************************

void out_ambig_table(const ISODUAL3D::AMBIG_TABLE & ambig_table)
{
  cout << "Dimension: " << DIM3 << endl;
  cout << "Number of table entries: "
       << ambig_table.NumTableEntries() << endl;
  for (int it = 0; it < ambig_table.NumTableEntries(); it++) {
    cout << "Entry " << it << ":";
    if (ambig_table.IsAmbiguous(it)) 
      { cout << "  Ambiguous."; }
    else
      { cout << "  Not ambiguous."; }
    cout << "  Num ambiguous facets: " 
         << int(ambig_table.NumAmbiguousFacets(it))
         << ".";
    if (ambig_table.NumAmbiguousFacets(it) > 0) {
      cout << "  Ambiguous facets:";
      for (int k = 0; k < NUM_CUBE_FACETS3D; k++) {
        if (ambig_table.IsFacetAmbiguous(it, k)) {
          cout << "  " << k;
        }
      }
    }
    cout << endl;

    cout << "    Num positive connected components: "
         << int(ambig_table.NumPosComponents(it))
         << "  Num negative connected components: "
         << int(ambig_table.NumNegComponents(it)) << "."
         << endl;
  }
  cout << endl;
}

// **************************************************
// Check routines
// **************************************************

void out_ambig_table_error
(const ISODUAL3D::AMBIG_TABLE & ambig_table,
 const int it, const bool flag_pos, const int expected_val)
{
  cerr << "Ambig table error, entry " << it << "." << endl;
  if (flag_pos) {
    cerr << "  Num pos components = " << ambig_table.NumPosComponents(0)
         << ".  Should be " << expected_val << "." << endl;
  }
  else {
    cerr << "  Num neg components = " << ambig_table.NumNegComponents(0)
         << ".  Should be " << expected_val << "." << endl;
  }
}

void out_ambig_table_zero_error
(const ISODUAL3D::AMBIG_TABLE & ambig_table,
 const int it, const bool flag_pos)
{
  cerr << "Ambig table error, entry " << it << "." << endl;
  if (flag_pos) {
    cerr << "  Num pos components = 0.  Should be non-zero." << endl;
  }
  else {
    cerr << "  Num neg components = 0.  Should be non-zero." << endl;
  }
}

void out_symmetry_error
(const ISODUAL3D::AMBIG_TABLE & ambig_table, const int it)
{
  int it2 = ambig_table.NumTableEntries() - it - 1;

  cerr << "Ambig table error, entry " << it 
       << " and complementary entry " << it2 << "."
       << endl;
  cerr << "  Num pos components[" << it << "] != Num neg components["
       << it2 << "]." << endl;
}

void check_ambig_table(const ISODUAL3D::AMBIG_TABLE & ambig_table)
{
  if (ambig_table.NumPosComponents(0) != 0) 
    { out_ambig_table_error(ambig_table, 0, true, 0); }
  if (ambig_table.NumNegComponents(0) != 1)
    { out_ambig_table_error(ambig_table, 0, false, 1); }
  if (ambig_table.NumPosComponents(255) != 1)
    { out_ambig_table_error(ambig_table, 255, true, 1); }
  if (ambig_table.NumNegComponents(255) != 0)
    { out_ambig_table_error(ambig_table, 255, false, 0); }

  for (int it = 1; it+1 < ambig_table.NumTableEntries(); it++) {
    if (ambig_table.NumPosComponents(it) == 0)
      { out_ambig_table_zero_error(ambig_table, it, true); }
    if (ambig_table.NumNegComponents(it) == 0)
      { out_ambig_table_zero_error(ambig_table, it, false); }
  }

  for (int it = 0; it < ambig_table.NumTableEntries(); it++) {
    int it2 = ambig_table.NumTableEntries() - it - 1;
    if (ambig_table.NumPosComponents(it) != 
        ambig_table.NumNegComponents(it2)) {
      out_symmetry_error(ambig_table, it);
    }

    if (ambig_table.NumPosComponents(it2) != 
        ambig_table.NumNegComponents(it)) {
      out_symmetry_error(ambig_table, it2);
    }
  }

}

// **************************************************
// Miscellaneous routines
// **************************************************

void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}

void usage_msg()
{
  cerr << "Usage: testambig <dimension>" << endl;
}

void usage_error()
{
  usage_msg();
  exit(10);
}

void parse_command_line(int argc, char **argv)
{

  int iarg = 1;
  while (iarg < argc && argv[iarg][0] == '-') {

    string s = string(argv[iarg]);
    usage_error();

    iarg++;
  }

  if (iarg != argc) { usage_error(); }
}

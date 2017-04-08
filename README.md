# LR-eigenvalues
Computes the eigenvalues of any kind of quadratic matrixes.

Usage
===
```
app.out [input_filename] [output_filename] [options]
```
Where _options_ include:
* `-d` print debug messages
* `-e` print errors
* `-p` print matrix
* `-t` print execution time
* `-prec=<num>` precision (default: 1e-14)
* `-eps=<num>` epsilon (a = 0 if a < eps; default: 1e-10)
* `-max_iter=<num>` maximum iterations (default: 0, i.e. no limit)
* `-h` or `-?` print help and exit

Default input filename is `infile.txt` and input filename is `outfile.txt`.

Input file format:
```
	n
   	a_1_1 a_1_2 ... a_1_n
   	a_2_1 a_2_2 ... a_2_n
   	.....................
   	a_n_1 a_n_2 ... a_n_n
```
Where `n` - input matrix size, `a_i_j` - matrix element on i-th row and j-th column.

Output file format (if there were errors, there will be just `0`):
```
	n
	e_1
	e_2
	...
	e_n
```
`e_i` - Eigenvalues, sorted ascending.

License
--------
	MIT License
	
	Copyright (c) 2017 exaltation
	
	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:
	
	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.
	
	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE.

/*
   Copyright 2015 Vinhthuy Phan
	Compressed FM index.
*/
package fmic

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io/ioutil"
	"os"
	"path"
	"strconv"
	"strings"
	"sync"
	"unsafe"
)

type Symb_OCC struct {
	Symb int
	OCC  []indexType
}

//-----------------------------------------------------------------------------
func check_for_error(e error) {
	if e != nil {
		panic(e)
	}
}

//-----------------------------------------------------------------------------
// Save the index to directory.

func _save_indexType(s []indexType, filename string) {
	f, err := os.Create(filename)
	check_for_error(err)
	defer f.Close()
	w := bufio.NewWriter(f)
	err = binary.Write(w, binary.LittleEndian, s)
	check_for_error(err)
	w.Flush()
}

func _save_sequenceType(s []sequenceType, filename string) {
	f, err := os.Create(filename)
	check_for_error(err)
	defer f.Close()
	w := bufio.NewWriter(f)
	err = binary.Write(w, binary.LittleEndian, s)
	check_for_error(err)
	w.Flush()
}

// ------------------------------------------------------------------
// save_option:
// 	0 - do not save suffix array and seq
//		1 - save suffix array, but not seq
//		2 - save both suffix array and seq
// ------------------------------------------------------------------
func (I *IndexC) SaveCompressedIndex(save_option int) {
	dir := I.input_file + ".fmi"
	os.Mkdir(dir, 0777)

	var wg sync.WaitGroup
	wg.Add(len(I.SYMBOLS) + 4)

	go func() {
		defer wg.Done()
		err := ioutil.WriteFile(path.Join(dir, "bwt"), I.BWT, 0666)
		check_for_error(err)
	}()

	go func() {
		defer wg.Done()
		_save_sequenceType(I.SSA, path.Join(dir, "ssa"))
	}()

	go func() {
		defer wg.Done()
		if save_option == 1 || save_option == 2 {
			_save_indexType(I.SA, path.Join(dir, "sa"))
		}
	}()

	go func() {
		defer wg.Done()
		if save_option == 2 {
			err := ioutil.WriteFile(path.Join(dir, "seq"), I.SEQ, 0666)
			check_for_error(err)
		}
	}()

	for symb := range I.OCC {
		go func(symb byte) {
			defer wg.Done()
			_save_indexType(I.OCC[symb], path.Join(dir, "occ."+string(symb)))
		}(symb)
	}

	f, err := os.Create(path.Join(dir, "others"))
	check_for_error(err)
	defer f.Close()
	w := bufio.NewWriter(f)
	fmt.Fprintf(w, "%d %d %d %d %t %d\n", I.LEN, I.OCC_SIZE, I.END_POS, I.M, I.Multiple, save_option)
	for i := 0; i < len(I.SYMBOLS); i++ {
		symb := byte(I.SYMBOLS[i])
		fmt.Fprintf(w, "%s %d %d %d\n", string(symb), I.Freq[symb], I.C[symb], I.EP[symb])
	}
	w.Flush()

	// save genome info
	f, err = os.Create(path.Join(dir, "genome_lengths"))
	check_for_error(err)
	defer f.Close()
	w = bufio.NewWriter(f)
	for i := 0; i < len(I.GENOME_ID); i++ {
		fmt.Fprintf(w, "%d %s %s\n", I.LENS[i], I.GENOME_ID[i], I.GENOME_DES[i])
	}
	w.Flush()

	wg.Wait()
}

// ------------------------------------------------------------------
// save_option:
// 	0 - suffix array and seq were not saved
//		1 - suffix array was saved; seq was not
//		2 - both suffix array and seq were saved
// ------------------------------------------------------------------
func LoadCompressedIndex(dir string) *IndexC {
	I := new(IndexC)

	// First, load "others"
	f, err := os.Open(path.Join(dir, "others"))
	check_for_error(err)
	defer f.Close()

	var symb byte
	var freq, c, ep indexType
	var save_option int
	scanner := bufio.NewScanner(f)
	scanner.Scan()
	fmt.Sscanf(scanner.Text(), "%d%d%d%d%t%d\n", &I.LEN, &I.OCC_SIZE, &I.END_POS, &I.M, &I.Multiple, &save_option)

	I.Freq = make(map[byte]indexType)
	I.C = make(map[byte]indexType)
	I.EP = make(map[byte]indexType)
	for scanner.Scan() {
		fmt.Sscanf(scanner.Text(), "%c%d%d%d", &symb, &freq, &c, &ep)
		I.SYMBOLS = append(I.SYMBOLS, int(symb))
		I.Freq[symb], I.C[symb], I.EP[symb] = freq, c, ep
	}

	// load genome_info
	f, err = os.Open(path.Join(dir, "genome_lengths"))
	check_for_error(err)
	defer f.Close()
	scanner = bufio.NewScanner(f)
	var items []string
	for scanner.Scan() {
		items = strings.SplitN(strings.TrimSpace(scanner.Text()), " ", 3)
		cur_len, _ := strconv.Atoi(items[0])
		I.GENOME_ID = append(I.GENOME_ID, items[1])
		I.GENOME_DES = append(I.GENOME_ID, items[2])
		I.LENS = append(I.LENS, indexType(cur_len))
	}

	// Second, load Suffix array, BWT and OCC
	I.OCC = make(map[byte][]indexType)
	var wg sync.WaitGroup
	wg.Add(len(I.SYMBOLS) + 4)

	go func() {
		defer wg.Done()
		I.BWT, err = ioutil.ReadFile(path.Join(dir, "bwt"))
		check_for_error(err)
	}()

	go func() {
		defer wg.Done()
		I.SSA = _load_sequenceType(path.Join(dir, "ssa"), I.LEN)
	}()

	go func() {
		defer wg.Done()
		if save_option == 1 || save_option == 2 {
			I.SA = _load_indexType(path.Join(dir, "sa"), I.LEN)
		}
	}()

	go func() {
		defer wg.Done()
		if save_option == 2 {
			I.SEQ, err = ioutil.ReadFile(path.Join(dir, "seq"))
			check_for_error(err)
		}
	}()

	Symb_OCC_chan := make(chan Symb_OCC)
	for _, symb := range I.SYMBOLS {
		go func(symb int) {
			defer wg.Done()
			Symb_OCC_chan <- Symb_OCC{symb, _load_indexType(path.Join(dir, "occ."+string(symb)), I.OCC_SIZE)}
		}(symb)
	}
	go func() {
		wg.Wait()
		close(Symb_OCC_chan)
	}()

	for symb_occ := range Symb_OCC_chan {
		I.OCC[byte(symb_occ.Symb)] = symb_occ.OCC
	}
	return I
}

//-----------------------------------------------------------------------------
func _load_indexType(filename string, length indexType) []indexType {
	f, err := os.Open(filename)
	check_for_error(err)
	defer f.Close()
	numBytes := uint(unsafe.Sizeof(indexType(0)))
	v := make([]indexType, length)

	scanner := bufio.NewScanner(f)
	scanner.Split(bufio.ScanBytes)
	for i, b := 0, uint(0); scanner.Scan(); b++ {
		if b == numBytes {
			b, i = 0, i+1
		}
		v[i] += indexType(scanner.Bytes()[0]) << (b * 8)
	}
	return v
}

//-----------------------------------------------------------------------------
func _load_sequenceType(filename string, length indexType) []sequenceType {
	f, err := os.Open(filename)
	check_for_error(err)
	defer f.Close()
	numBytes := uint(unsafe.Sizeof(sequenceType(0)))
	v := make([]sequenceType, length)

	scanner := bufio.NewScanner(f)
	scanner.Split(bufio.ScanBytes)
	for i, b := 0, uint(0); scanner.Scan(); b++ {
		if b == numBytes {
			b, i = 0, i+1
		}
		v[i] += sequenceType(scanner.Bytes()[0]) << (b * 8)
	}
	return v
}

//-----------------------------------------------------------------------------

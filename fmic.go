/*
   Copyright 2015 Vinhthuy Phan
	Reversed Compressed FM index.
*/
package fmic

import (
	"bufio"
	"bytes"
	"fmt"
	"math"
	"math/rand"
	"os"
	"sort"
)

//-----------------------------------------------------------------------------
// Global variables: sequence (SEQ), suffix array (SA), BWT, FM index (C, OCC)
//-----------------------------------------------------------------------------

type IndexC struct {
	SEQ []byte
	BWT []byte
	SA  []indexType          // suffix array
	SSA []sequenceType       // SSA[i] stores the sequence containing position SA[i]
	C   map[byte]indexType   // count table
	OCC map[byte][]indexType // occurence table

	END_POS indexType          // position of "$" in the text
	SYMBOLS []int              // sorted symbols
	EP      map[byte]indexType // ending row/position of each symbol

	LEN        indexType
	LENS       []indexType
	GENOME_ID  []string
	GENOME_DES  []string
	OCC_SIZE   indexType
	Freq       map[byte]indexType // Frequency of each symbol
	M          int                // Compression ratio
	Multiple   bool               // True if the input contains multiple sequences
	input_file string
}

//-----------------------------------------------------------------------------
// Build FM index given the file storing the text.
// multiple is true if the input file contains multiple sequences
// compression ratio >=1
//-----------------------------------------------------------------------------
func CompressedIndex(file string, multiple bool, compression_ratio int) *IndexC {
	I := new(IndexC)
	I.input_file = file
	I.M = compression_ratio
	I.Multiple = multiple

	// GET THE SEQUENCE
	I.ReadFasta(file)

	// BUILD SUFFIX ARRAY
	I.LEN = indexType(len(I.SEQ))
	I.OCC_SIZE = indexType(math.Ceil(float64(I.LEN/indexType(I.M)))) + 1
	I.SA = make([]indexType, I.LEN)
	var SID []sequenceType
	if I.Multiple {
		I.SSA = make([]sequenceType, I.LEN)
		SID = make([]sequenceType, I.LEN)
	}
	SA := make([]int, I.LEN)
	ws := &WorkSpace{}
	ws.ComputeSuffixArray(I.SEQ, SA)
	sid := sequenceType(0)
	for i := range SA {
		I.SA[i] = indexType(SA[i])
		if I.Multiple {
			SID[i] = sid
			if I.SEQ[i] == '|' {
				sid++
			}
		}
	}

	// BUILD BWT
	I.Freq = make(map[byte]indexType)
	I.BWT = make([]byte, I.LEN)
	var i indexType
	for i = 0; i < I.LEN; i++ {
		I.Freq[I.SEQ[i]]++
		if I.SA[i] == 0 {
			I.BWT[i] = I.SEQ[I.LEN-1]
		} else {
			I.BWT[i] = I.SEQ[I.SA[i]-1]
		}
		if I.BWT[i] == '$' {
			I.END_POS = i
		}
		if I.Multiple {
			// I.SSA[i] = SID[I.SA[i]]
			I.SSA[i] = sid - SID[I.SA[i]]   // This is because I.SEQ is reversed.
		}
	}

	// BUILD COUNT AND OCCURENCE TABLE
	I.C = make(map[byte]indexType)
	I.OCC = make(map[byte][]indexType)
	for c := range I.Freq {
		I.SYMBOLS = append(I.SYMBOLS, int(c))
		I.OCC[c] = make([]indexType, I.OCC_SIZE)
		I.C[c] = 0
	}
	sort.Ints(I.SYMBOLS)
	I.EP = make(map[byte]indexType)
	count := make(map[byte]indexType)

	for j := 1; j < len(I.SYMBOLS); j++ {
		curr_c, prev_c := byte(I.SYMBOLS[j]), byte(I.SYMBOLS[j-1])
		I.C[curr_c] = I.C[prev_c] + I.Freq[prev_c]
		I.EP[curr_c] = I.C[curr_c] + I.Freq[curr_c] - 1
		count[curr_c] = 0
	}

	for j := 0; j < len(I.BWT); j++ {
		count[I.BWT[j]] += 1
		if j%I.M == 0 {
			for symbol := range I.OCC {
				I.OCC[symbol][int(j/I.M)] = count[symbol]
			}
		}
	}

	return I
}

//-----------------------------------------------------------------------------
func (I *IndexC) Occurence(c byte, pos indexType) indexType {
	i := indexType(pos / indexType(I.M))
	count := I.OCC[c][i]
	for j := i*indexType(I.M) + 1; j <= pos; j++ {
		if I.BWT[j] == c {
			count += 1
		}
	}
	return count
}

// -----------------------------------------------------------------------------
// Returns starting, ending positions (sp, ep) and last-matched position (i)

func (I *IndexC) Search(query []byte) (int, int) {
	var offset indexType
	var i int
	start_pos := 0
	c := query[start_pos]
	sp, ok := I.C[c]
	if !ok {
		panic("Unknown character: " + string(c))
	}
	ep := I.EP[c]
	// fmt.Println(i, string(c), sp, ep)
	for i = int(start_pos + 1); sp <= ep && i < len(query); i++ {
		c = query[i]
		offset, ok = I.C[c]
		if !ok {
			panic("Unknown character: " + string(c))
		}
		sp = offset + I.Occurence(c, sp-1)
		ep = offset + I.Occurence(c, ep) - 1
		// fmt.Println(i, string(c), sp, ep)
	}
	return int(sp), int(ep)
}

//-----------------------------------------------------------------------------
func (I *IndexC) TestSearch(query []byte) map[sequenceType][]indexType {
	const maxSize = 1
	idSet := map[sequenceType][]indexType{}
	var offset indexType
	var i int
	for start_pos := 0; start_pos < len(query); start_pos = i {
		flag := true
		c := query[start_pos]
		sp, ok := I.C[c]
		ep := I.EP[c]
		for i = start_pos; sp <= ep && i < len(query); i++ {
			if ep-sp <= maxSize && flag == true {
				flag = false
				for i := sp; i <= ep; i++ {
					idSet[I.SSA[i]] = append(idSet[I.SSA[i]], I.SA[i])
				}
			}
			c = query[i]
			offset, ok = I.C[c]
			if !ok {
				return idSet
			}
			sp = offset + I.Occurence(c, sp-1)
			ep = offset + I.Occurence(c, ep) - 1
		}
	}
	return idSet
}

//-----------------------------------------------------------------------------
func (I *IndexC) regionSearch(query []byte, start_pos int) (int, int, map[sequenceType]indexType) {
	const maxSize = 10
	idSet := map[sequenceType]indexType{}
	flag := true
	var id sequenceType
	var offset, pos indexType
	var i int
	c := query[start_pos]
	sp, ok := I.C[c]
	if !I.Multiple || !ok {
		return -1,-1,idSet
	}
	ep := I.EP[c]
	for i = int(start_pos + 1); sp < ep && i < len(query); i++ {
		if ep-sp <= maxSize && flag == true {
			flag = false
			for i := sp; i <= ep; i++ {
				idSet[I.SSA[i]] = I.SA[i]
			}
			// If all regions are the same, return.  Else, continue.
			if len(idSet) == 1 {
				for id, pos = range idSet {
					return int(id), int(pos), idSet
				}
			}
		}
		c = query[i]
		offset, ok = I.C[c]
		if !ok {
			return -1,-1,idSet
		}
		sp = offset + I.Occurence(c, sp-1)
		ep = offset + I.Occurence(c, ep) - 1
	}
	if sp == ep {
		idSet[I.SSA[sp]] = I.SA[sp]
		return int(I.SSA[sp]), int(I.SA[sp]), idSet
	} else {
		return -1,-1,idSet
	}
}

//-----------------------------------------------------------------------------
func (I *IndexC) FindGenomeD(query1 []byte, query2 []byte, maxInsert int) map[int]int {
	id1, pos1, idSet1 := I.regionSearch(query1, 0)
	id2, pos2, idSet2 := I.regionSearch(query2, 0)
	out := map[int]int{}
	// fmt.Println("\t",id1,pos1,idSet1,"\t",id2,pos2,idSet2)
	if id1==id2 && id1!=-1 && ((pos1>=pos2 && int(pos1-pos2)<=maxInsert)||(pos2>pos1 && int(pos2-pos1)<=maxInsert)) {
		out[int(id1)] = 1
	} else {
		for id, p1 := range idSet1 {
			if p2, ok := idSet2[id]; ok && int(id)!=-1 {
				if (p1 >= p2 && int(p1-p2) <= maxInsert) || (p2 > p1 && int(p2-p1) <= maxInsert) {
					out[int(id)] = 1
				}
			}
		}
	}
	return out
}

//-----------------------------------------------------------------------------
func (I *IndexC) FindGenomeR(query1 []byte, query2 []byte, maxInsert int, rounds int) map[int]int {
	k1, k2 := 0,0	// init round starts from fixed index
	end := 20
	regions := map[int]int{}
	for i:=0; i<rounds; i++ {
		id1, pos1, idSet1 := I.regionSearch(query1, k1)
		id2, pos2, idSet2 := I.regionSearch(query2, k2)
		out := map[int]int{}
		if id1==id2 && id1!=-1 && ((pos1>=pos2 && int(pos1-pos2)<=maxInsert)||(pos2>pos1 && int(pos2-pos1)<=maxInsert)) {
			// fmt.Println("1:", pos1, pos2, out)
			out[int(id1)] = 1
			return out
		} else {
			for id, p1 := range idSet1 {
				if p2, ok := idSet2[id]; ok && int(id)!=-1 {
					if (p1 >= p2 && int(p1-p2) <= maxInsert) || (p2 > p1 && int(p2-p1) <= maxInsert) {
						out[int(id)] = 1
						regions[int(id)]++
					}
				}
			}
			if len(out) == 1 {  // conservative
				// fmt.Println("2:", pos1, pos2, out)
				return out
			}
		}
		k1 = rand.Intn(len(query1)-end)
		k2 = rand.Intn(len(query2)-end)
	}
	// fail
	// reg, max := -1, 0
	// for r,count := range regions {
	// 	if count > max {
	// 		reg, max = r, count
	// 	}
	// }
	// fmt.Println(">>>>Highest count is", reg, max)
	// fmt.Println(">>>", regions)
	return map[int]int{}
}
//-----------------------------------------------------------------------------
// func (I *IndexC) FindGenome(query1 []byte, query2 []byte, randomized_round, maxInsert int) map[int]int {
// 	var gid1, gid2 map[sequenceType]indexType
// 	var pos int
// 	for i := 0; i < randomized_round; i++ {
// 		pos = rand.Intn(len(query1)-10)
// 		gid1 = I.regionSearch(query1, pos)
// 		pos = rand.Intn(len(query2)-10)
// 		gid2 = I.regionSearch(query2, pos)
// 		out := make(map[int]int)
// 		for gid, p1 := range gid1 {
// 			if p2, ok := gid2[gid]; ok {
// 				if (p1 >= p2 && int(p1-p2) <= maxInsert) || (p2 > p1 && int(p2-p1) <= maxInsert) {
// 					out[int(gid)] = 1
// 				}
// 			}
// 		}
// 		if len(out) > 0 {
// 			// fmt.Println("\t", len(gid1), len(gid2), len(out))
// 			return out
// 		}
// 	}
// 	return map[int]int{}
// }

//-----------------------------------------------------------------------------
// Guess which sequence contains the query.
// If randomized_round is 0, there is no randomization. The search begins at the.
//-----------------------------------------------------------------------------

// func (I *IndexC) Guess(query []byte, randomized_round int) (int, int) {
// 	var seq, count int
// 	// var start_pos, end_pos int
// 	if randomized_round == 0 {
// 		seq, count, _ = I._guess(query, len(query)-1)
// 		return seq, count
// 	} else {
// 		for i := 0; i < randomized_round; i++ {
// 			// start_pos = rand.Intn(len(query))
// 			// seq, count, end_pos = I._guess(query, start_pos)
// 			// fmt.Println(end_pos, start_pos, "<")
// 			seq, count, _ = I._guess(query, rand.Intn(len(query)))
// 			if seq >= 0 {
// 				return seq, count
// 			}
// 		}
// 		return -1, 0
// 	}
// }

// //-----------------------------------------------------------------------------
// func (I *IndexC) GuessPairD(query1 []byte, query2 []byte) int {
// 	var seq1, seq2, p1, p2 int
// 	max := len(query1)
// 	if max > len(query2) {
// 		max = len(query2)
// 	}
// 	maxInsert := 1500
// 	for pos := 15; pos < max; pos++ {
// 		seq1, _, p1 = I._guess(query1, pos)
// 		seq2, _, p2 = I._guess(query2, pos)

// 		// fmt.Println(seq1, p1, int(I.LEN)-p1+1, "|", seq2, p2, int(I.LEN)-p2+1)
// 		if seq1 == seq2 && seq1 != -1 &&
// 			((p1 >= p2 && p1-p2 <= maxInsert) || (p2 > p1 && p2-p1 <= maxInsert)) {
// 			return seq1
// 		}
// 	}
// 	return -1
// }

// //-----------------------------------------------------------------------------
// func (I *IndexC) GuessPair(query1 []byte, query2 []byte, randomized_round, maxInsert int) int {
// 	var seq1, seq2, p1, p2, pos int
// 	// var c1, c2 int
// 	for i := 0; i < randomized_round; i++ {
// 		pos = 10 + rand.Intn(len(query1)-10)
// 		// fmt.Printf("left ")
// 		seq1, _, p1 = I._guess(query1, pos)
// 		pos = 10 + rand.Intn(len(query2)-10)
// 		// fmt.Printf("right ")
// 		seq2, _, p2 = I._guess(query2, pos)

// 		// fmt.Println(seq1, p1, int(I.LEN)-p1+1, "|", seq2, p2, int(I.LEN)-p2+1)
// 		// fmt.Println(seq1, seq2, "\t", c1, c2, "\t", p1, p2)
// 		if seq1 == seq2 && seq1 != -1 &&
// 			((p1 >= p2 && p1-p2 <= maxInsert) || (p2 > p1 && p2-p1 <= maxInsert)) {
// 			return seq1
// 		}
// 	}
// 	return -1
// }

// //-----------------------------------------------------------------------------
// func (I *IndexC) _guess(query []byte, start_pos int) (int, int, int) {
// 	if !I.Multiple {
// 		return 0, -1, -1
// 	}
// 	var offset indexType
// 	var i int
// 	c := query[start_pos]
// 	sp, ok := I.C[c]
// 	if !ok {
// 		return -2, 0, 0
// 	}
// 	ep := I.EP[c]
// 	// fmt.Println(ep-sp+1, "\t", i, string(c), len(query))
// 	for i = int(start_pos - 1); sp < ep && i >= 0; i-- {
// 		c = query[i]
// 		offset, ok = I.C[c]
// 		if !ok {
// 			return -2, 0, 0
// 		}
// 		// if ep-sp <= 10 {
// 		// 	for k := sp; k <= ep; k++ {
// 		// 		fmt.Printf("%d ", I.SSA[k])
// 		// 	}
// 		// 	fmt.Println("[", start_pos-i, "]")
// 		// }
// 		sp = offset + I.Occurence(c, sp-1)
// 		ep = offset + I.Occurence(c, ep) - 1
// 		// fmt.Println(ep-sp+1, "\t", i, string(c), len(query))
// 	}
// 	if sp <= ep {
// 		for j := sp + 1; j <= ep; j++ {
// 			if I.SSA[j] != I.SSA[sp] {
// 				return -1, int(ep - sp + 1), -1
// 			}
// 		}
// 		return int(I.SSA[sp]), int(ep - sp + 1), int(I.SA[sp])
// 	} else {
// 		return -1, int(ep - sp + 1), -1
// 	}
// }

//-----------------------------------------------------------------------------
func (I *IndexC) ReadFasta(file string) {
	f, err := os.Open(file)
	check_for_error(err)
	defer f.Close()

	if file[len(file)-6:] != ".fasta" {
		panic("ReadFasta:" + file + "is not a fasta file.")
	}

	scanner := bufio.NewScanner(f)
	byte_array := make([]byte, 0)
	i := 0
	cur_len := 0
	for scanner.Scan() {
		line := scanner.Bytes()
		if len(line) > 0 {
			line = bytes.Trim(line, "\n\r ")
			if line[0] != '>' {
				byte_array = append(byte_array,line...)
				cur_len += len(line)
			} else {
				items := bytes.SplitN(line[1:], []byte{' '}, 2)
				I.GENOME_ID = append(I.GENOME_ID, string(items[0]))
				I.GENOME_DES = append(I.GENOME_DES, string(items[1]))
				if cur_len != 0 {
					I.LENS = append(I.LENS, indexType(cur_len))
				}
				cur_len = 0
				if len(byte_array) > 0 {
					byte_array = append(byte_array, byte('|'))
				}
			}
			i++
		}
	}
	I.LENS = append(I.LENS, indexType(cur_len))
	// Reverse the sequence
	for left, right := 0, len(byte_array)-1; left < right; left, right = left+1, right-1 {
	    byte_array[left], byte_array[right] = byte_array[right], byte_array[left]
	}
	I.SEQ = append(byte_array, byte('$'))
}

//-----------------------------------------------------------------------------
func (I *IndexC) Show() {
	fmt.Printf(" %6s %6s  OCC\n", "Freq", "C")
	for i := 0; i < len(I.SYMBOLS); i++ {
		c := byte(I.SYMBOLS[i])
		fmt.Printf("%c%6d %6d  %d\n", c, I.Freq[c], I.C[c], I.OCC[c])
	}
	fmt.Printf("SA ")
	for i := 0; i < len(I.SA); i++ {
		fmt.Print(I.SA[i], " ")
	}
	fmt.Printf("\nBWT ")
	for i := 0; i < len(I.BWT); i++ {
		fmt.Print(string(I.BWT[i]))
	}
	fmt.Println()
	fmt.Printf("\nSSA ")
	for i := 0; i < len(I.SSA); i++ {
		fmt.Printf("%d ", I.SSA[i])
	}
	fmt.Println()
	fmt.Println("SEQ", string(I.SEQ))
	for i:=0; i<len(I.SA); i++ {
		fmt.Printf("%4d %s\n", i, string(I.SEQ[I.SA[i]:]))
	}
}

//-----------------------------------------------------------------------------
func (I *IndexC) Check() {
	fmt.Println("Checking...")
	for i := 0; i < len(I.SYMBOLS); i++ {
		c := byte(I.SYMBOLS[i])
		fmt.Printf("%c%6d %6d  [", c, I.Freq[c], I.C[c])
		for j := 0; j < int(I.LEN); j++ {
			fmt.Printf("%d ", I.Occurence(c, indexType(j)))
		}
		fmt.Printf("]\n")
	}
	if len(I.SEQ) > 0 {
		S := make([]byte, len(I.SEQ))
		for i:=0; i<len(S); i++ {
			S[i] = I.SEQ[len(S) - i - 1]
		}
		S = S[1:]
		fmt.Println(string(S))
		ep, sp := I.Search(S)
		fmt.Println("Search for SEQ returns", sp, ep)
	}
}

//-----------------------------------------------------------------------------

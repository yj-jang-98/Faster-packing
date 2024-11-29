package utils

import (
	"bufio"
	"encoding/csv"
	"os"
)

func DataExport(data [][]float64, fileName string) {
	file, err := os.Create(fileName)
	if err != nil {
		panic(err)
	}
	wr := csv.NewWriter(bufio.NewWriter(file))
	wr.WriteAll(MatToString(data))
}

PDF = $(wildcard *.pdf)
PNG = $(PDF:%.pdf=%.png)

all: $(PNG)

%.png: %.pdf
	convert -flatten -density 500 $< -resize 25% $@

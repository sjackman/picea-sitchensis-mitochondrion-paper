pandoc_opt=-Fpandoc-crossref -Fpandoc-citeproc

.DELETE_ON_ERROR:
.SECONDARY:

all: psitchensismt.pdf

clean:
	rm -f psitchensismt.html psitchensismt.pdf psitchensismt-supp.html psitchensismt-supp.pdf

# Download the citation style language (CSL)
psitchensismt.csl:
	curl -o $@ https://www.zotero.org/styles/genome-biology-and-evolution

# Render Markdown to HTML using Pandoc
%.html: %.md
	pandoc $(pandoc_opt) -s -o $@ $<

# Render Markdown to PDF using Pandoc
%.pdf: %.md
	pandoc $(pandoc_opt) --pdf-engine=xelatex -o $@ $<

# Generate Table of Contents for supplemental material only
psitchensismt-supp.pdf: psitchensismt-supp.md
	pandoc $(pandoc_opt) --toc -o $@ $<

# Render Markdown to DOCX using Pandoc
%.docx: %.md
	pandoc $(pandoc_opt) -o $@ $<

# Render RMarkdown to HTML using R
%.html: %.rmd
	RScript -e 'rmarkdown::render("$<")'

# Fetch BibTex records from a list of DOI.
%.doi.bib: %.doi
	brew cite $$(<$<) | sort >$@

# Concatentate the citations with and without DOI.
%.bib: %.doi.bib %.other.bib
	sort $^ | sed 's~http://dx.doi.org~https://doi.org~' >$@

psitchensismt.docx: psitchensismt.bib psitchensismt.csl
psitchensismt.html: psitchensismt.bib psitchensismt.csl
psitchensismt.pdf: psitchensismt.bib psitchensismt.csl

psitchensismt-supp.docx: psitchensismt.bib psitchensismt.csl
psitchensismt-supp.html: psitchensismt.bib psitchensismt.csl
psitchensismt-supp.pdf: psitchensismt.bib psitchensismt.csl

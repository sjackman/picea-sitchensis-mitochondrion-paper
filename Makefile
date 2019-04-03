pandoc_opt=-Fpandoc-crossref -Fpandoc-citeproc

.DELETE_ON_ERROR:
.SECONDARY:

all: psitchensismt.pdf

clean:
	rm -f psitchensismt.html psitchensismt.pdf

# Download the citation style language (CSL)
psitchensismt.csl:
	curl -o $@ https://www.zotero.org/styles/genome-biology-and-evolution

# Render Markdown to HTML using Pandoc
%.html: %.md
	pandoc $(pandoc_opt) -s -o $@ $<

# Render Markdown to PDF using Pandoc
%.pdf: %.md
	pandoc $(pandoc_opt) --pdf-engine=xelatex -o $@ $<

# Render Markdown to DOCX using Pandoc
%.docx: %.md
	pandoc $(pandoc_opt) -o $@ $<

# Fetch BibTex records from a list of DOI.
%.doi.bib: %.doi
	brew cite $$(<$<) | sort >$@

# Concatentate the citations with and without DOI.
# Preserve title case.
%.bib: %.doi.bib %.other.bib
	sort $^ \
	| sed -E \
		-e 's/title={([^}]*)},/title={{\1}},/' \
		-e 's~http://dx.doi.org~https://doi.org~' \
		>$@

psitchensismt.docx: psitchensismt.bib psitchensismt.csl
psitchensismt.html: psitchensismt.bib psitchensismt.csl
psitchensismt.pdf: psitchensismt.bib psitchensismt.csl

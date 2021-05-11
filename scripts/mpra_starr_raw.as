table mpra_starr_raw
"BED6+2 MPRA/STARR-seq common file format for raw input (DNA) or output (RNA) count data. Each row should report on the unit closest to individual molecule - raw region, region-UMI, region-reaction combination etc. rawReads column would therefore ideally represent PCR duplicates (e.g. number of times the same UMI was sequenced). Names in column 4 should ideally NOT be unique and should represent unique DNA or RNA sequences (CRS, oligos, elements, fragments...)"
(
string  chrom;		"Reference sequence chromosome or scaffold"
uint    chromStart;	"Start position in chromosome"
uint    chromEnd;	"End position in chromosome"
string  name;		"Name of tested sequence or region"
uint    score;		"Set to 0, ignored"
char[1] strand;		"+ or - for strand, . for unknown"
uint    rawReads;	"Raw number of DNA or RNA reads containing this barcode or UMI (MPRA, eSTARR), or mapping to this pre-binned region, or region-UMI/reaction combination (STARR)"
string  barcode;	"barcode, UMI or reaction name (if used, else enter .)"
)
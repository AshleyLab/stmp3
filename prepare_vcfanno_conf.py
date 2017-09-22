#written by Noah Friedman
#writes a .toml configuration file for vcf anno based on inputs

confHead = '[[annotation]]\n'

confFooter = '[[postannotation]]\nfields=[\"lua_start\"]\nop=\"lua:lua_start - 2\"\nname=\"lua_start_minus_2\"\ntype=\"Integer\"\n'

#pass the function the open file and it'll write the appropriate line
def write_file_name(f, filename):
	line = 'file=' + '\"' + filename + '\"' + '\n'
	print line
	f.write(line) 

def write_annotation_fields(f, fields):
	writeStr = 'fields = ['
	cntr = 1
	for field in fields:
		s = '\"' + field + '\"'
		writeStr += s
		if cntr != len(fields):
			writeStr += ', '
	writeStr += ']\n'
	print writeStr
	f.write(writeStr)

def write_ops_line(f, fields):
	writeStr = 'ops=['
	writeStr += ','.join(['\"first\"' for i in range(len(fields))])
	writeStr += ']\n'
	print writeStr
	f.write(writeStr)

def write_conf_file(fname, annotationFileDict):
	with open(fname, 'w') as f:
		for key, value in annotationFileDict.items():
			f.write(confHead)
			write_file_name(f, key)
			write_annotation_fields(f, value)
			write_ops_line(f, value)
			f.write('\n')
		f.write(confFooter)











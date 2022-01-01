extern "C" {
#include "frustrum_ifails.h"
#include "frustrum_tokens.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
}

#ifdef _WIN32
extern "C" {
#include <Windows.h>
}
#endif

#ifndef _WIN32
#define memcpy_s( dest, dest_size, src, len) memcpy(dest, src, len)
#define strcpy_s( dest, dest_size, src)      strcpy(dest, src)
#define strcat_s( dest, dest_size, src)      strcat(dest, src)
#endif

#include <map>
#include <fstream>
#include <string>
#include <memory>

static int frustrum_starts = 0;
static int num_open_files = 0;

std::string end_of_header_s = "**END_OF_HEADER";
static char end_of_header[32] = "**END_OF_HEADER";

static char newline_c = '\n';

class PSFile {
public:
	PSFile(const char* name, const int* namelen) {
		if (*namelen > 2 && name[0] == '*' && name[1] == '*') {
			data.append(name);
			loc = 0;
			is_text = true;
		}
		else {
			file.open(name);
			is_text = false;
		}
	}

	void skip_header() {
		if (is_text) {
			skip_header_text();
		}
		else {
			skip_header_file();
		}
	}

	void read(int max, char* buffer, int* n_read, int* ifail) {
		if (is_text) {
			read_text(max, buffer, n_read, ifail);
		}
		else {
			read_file(max, buffer, n_read, ifail);
		}
	}

private:

	void skip_header_text() {
		size_t end_of_header_start = data.find(end_of_header);
		loc = data.find("\n", end_of_header_start);
	}

	void skip_header_file() {
		std::string line;
		while (file.good() && line.compare(0, end_of_header_s.length(), end_of_header_s)) {
			std::getline(file, line);
		}
		if (!file.good()) {
			file.clear();
			file.seekg(0);
		}
	}

	void read_text(int max, char* buffer, int* n_read, int* ifail) {
		*ifail = FR_no_errors;
		if (loc >= data.size()) {
			*ifail = FR_end_of_file;
		}
		if (loc + max > data.size()) {
			*n_read = data.size() - loc;
		}
		else {
			*n_read = max;
		}
		std::string to_return = data.substr(loc, *n_read);
		loc += *n_read;
		memcpy_s(buffer, max, to_return.c_str(), *n_read);
	}

	void read_file(int max, char* buffer, int* n_read, int* ifail) {
		/*
		// No idea why this failed -- it was producing garbage
		// after a certain point in the file.
		if (file.good()) {
			*n_read = (int)file.readsome(buffer, max);
			*ifail = FR_no_errors;
		}
		else {
			*ifail = FR_end_of_file;
		}
		*/
		*ifail = FR_no_errors;
		*n_read = 0;
		for (int i = 0; i < max; ++i) {
			if (file.good()) {
				buffer[i] = file.get();
				++(*n_read);
				if (buffer[i] == newline_c) break;
			}
			else {
				if (*n_read < max) {
					*ifail = FR_end_of_file;
				}
				break;
			}
		}
	}

	bool is_text;
	std::ifstream file;
	std::string data;
	size_t loc;
};

std::map<int, std::shared_ptr<PSFile>> open_files;
int next_file_id = 0;

extern void FSTART(int* ifail)
{
	*ifail = FR_unspecified;
	if (frustrum_starts == 0) {
		// Setup any data structures
	}
	++frustrum_starts;
	*ifail = FR_no_errors;
}

extern void FSTOP(int* ifail)
{
	*ifail = FR_unspecified;
	if (frustrum_starts <= 0) return;
	--frustrum_starts;
	if (frustrum_starts == 0) {
		next_file_id = 0;
		open_files.clear();
	}
	*ifail = FR_no_errors;
}

extern void FMALLO(int* nbytes, char** memory, int* ifail)
{
	*ifail = FR_unspecified;
	if (frustrum_starts <= 0) {
		*memory = 0;
		return;
	}
	*memory = (char*)malloc(*nbytes);
	if (*memory == NULL) {
		*ifail = FR_memory_full;
		return;
	}
	*ifail = FR_no_errors;
}

extern void FMFREE(int* nbytes, char** memory, int* ifail)
{
	*ifail = FR_unspecified;
	if (frustrum_starts <= 0) {
		return;
	}
	free(*memory);
	*ifail = FR_no_errors;
}

extern void FFOPRD(const int* guise, const int* format, const char* name,
	const int* namlen, const int* skiphd, int* strid, int* ifail)
{
	*ifail = FR_unspecified;
	*strid = -1;
	if (frustrum_starts <= 0) return;
	auto file = std::make_shared<PSFile>(name, namlen);
	open_files.emplace(next_file_id, file);
	if (*skiphd) {
		open_files[next_file_id]->skip_header();
	}
	*strid = next_file_id;
	++next_file_id;
	*ifail = FR_no_errors;
}

extern void FFOPWR(const int* guise, const int* format, const char* name,
	const int* namlen, const char* pd2hdr, const int* pd2len,
	int* strid, int* ifail)
{
	// Dummy Function - we don't ever write
	*strid = 1;
	*ifail = FR_no_errors;
}

extern void FFWRIT(const int* guise, const int* strid, const int* nchars,
	const char* bufer, int* ifail)
{
	// Dummy Function - we don't ever write
	*ifail = FR_no_errors;
}

extern void FFREAD(const int* guise, const int* strid, const int* nmax,
	char* buffer, int* nactual, int* ifail)
{
	*ifail = FR_unspecified;
	*nactual = 0;
	if (frustrum_starts <= 0) return;
	auto file = open_files.find(*strid);
	if (file != open_files.end()) {
		*ifail = FR_no_errors;
		file->second->read(*nmax, buffer, nactual, ifail);
	}
}

extern void FFCLOS(const int* guise, const int* strid, const int* action,
	int* ifail)
{
	*ifail = FR_unspecified;
	if (frustrum_starts <= 0) return;
	auto it = open_files.find(*strid);
	*ifail = FR_no_errors;
	if (open_files.erase(*strid) == 0) {
		*ifail = FR_close_fail;
	}
}
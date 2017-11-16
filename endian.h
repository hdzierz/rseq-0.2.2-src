#ifndef ENDIAN_H
#define ENDIAN_H

#include "system_utils.h"

static inline void reverse_buf(char *buf, const int size){
	int i;
	for (i = 0; i < size / 2; i++) {
		char temp = buf[i];
		buf[i] = buf[size - 1 - i];
		buf[size - 1 - i] = temp;
	}
}

inline size_t big_endian_fread( void * ptr, size_t size, size_t count, FILE * stream ){
	char *temp;
	temp = new char[size * count];
	size_t result = fread(temp, size, count, stream);
	if (is_little_endian()) {
		int i;
		for (i = 0; i < (int)count; i++) {
			reverse_buf(temp + (int)size * i, (int)size);
		}
	}
	memcpy(ptr, temp, size*count);
	delete[] temp;
	return result;
}

inline size_t big_endian_fwrite( const void * ptr, size_t size, size_t count, FILE * stream ){
	char *temp;
	temp = new char[size * count];
	memcpy(temp, ptr, size*count);
	if (is_little_endian()) {
		int i;
		for (i = 0; i < (int)count; i++) {
			reverse_buf(temp + (int)size * i, (int)size);
		}
	}
	size_t result = fwrite(temp, size, count, stream);
	delete[] temp;
	return result;
}

inline size_t little_endian_fread( void * ptr, size_t size, size_t count, FILE * stream ){
	char *temp;
	temp = new char[size * count];
	size_t result = fread(temp, size, count, stream);
	if (!is_little_endian()) {
		int i;
		for (i = 0; i < (int)count; i++) {
			reverse_buf(temp + (int)size * i, (int)size);
		}
	}
	memcpy(ptr, temp, size*count);
	delete[] temp;
	return result;
}

inline size_t little_endian_fwrite( const void * ptr, size_t size, size_t count, FILE * stream ){
	char *temp;
	temp = new char[size * count];
	memcpy(temp, ptr, size*count);
	if (!is_little_endian()) {
		int i;
		for (i = 0; i < (int)count; i++) {
			reverse_buf(temp + (int)size * i, (int)size);
		}
	}
	size_t result = fwrite(temp, size, count, stream);
	delete[] temp;
	return result;
}

inline void affy_bug_fix(float &data){
	reverse_buf((char*)&data, 4);
	int temp = (int)data;
	reverse_buf((char*)&temp, 4);
	data = (float)temp;
}

#endif //ENDIAN_H

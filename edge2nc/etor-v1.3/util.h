/*----------------------------------------------------------------------**
**
** util.h - Header file for util library
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
**
** USAGE:
**
** PROCESSING:
**
**
** COPYRIGHT NOTICE
**
**	Copyright (c) 1997 by Enterprise Electronics Corporation
**	All Rights Reserved
** 
** This program is  copyright  by  Enterprise  Electronics  Corpora-
** tion,    Enterprise,  Alabama,  USA  36330 (334) 347-3478.  It is
** licensed for  use  on  a  specific  CPU   and   is  not    to  be
** transferred  or otherwise divulged.   Copies  or modifications of
** this program must carry this copyright notice.
** 
**
**
** HEADER INFOMATION
**
**	Software Suite 		- EDGE
**	Package			- UTIL
**	Reference number	- SP1/PGM/
**	Revision number		- $Revision: 1.1 $
**	Release State		- $State: Exp $
**	Author, designer	- Bob Stafford
** Modification Date		- $Date: 1999/03/31 22:35:30 $
** Modified by			- $Author: merritt $
** $Source: /nfs/trmm/src/CVS/etor/util.h,v $
**
** MODIFICATION RECORD
**
** $Log: util.h,v $
** Revision 1.1  1999/03/31 22:35:30  merritt
** the etor lib
**
**
**----------------------------------------------------------------------*/

/*
** Include the prototype file
*/

#define GETENVINT(E,V) { char *temp=getenv(E); if (temp) *(V)=atoi(temp);}
#define GETENVFLOAT(E,V) { char *temp=getenv(E); if (temp) *(V)=atof(temp);}
#define GETENVSTRING(E,V) { char *temp=getenv(E); if (temp) strcpy((V),temp);}
#define GETENVSTRINGP(E,V) { char *temp=getenv(E); if (temp) *(V)=temp;}

/*-END OF MODULE--------------------------------------------------------*/

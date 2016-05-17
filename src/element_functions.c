#include <stdlib.h>
#include <stdio.h>

#include "database.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Simple element-related functions:
 *			int        is_ELEMENT_present(const unsigned int type);
 *			*S_ELEMENT get_ELEMENT_type(const unsigned int type);
 *			*S_ELEMENT get_ELEMENT_Eclass(const unsigned int Eclass, const unsigned int Esubclass);
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

unsigned int is_ELEMENT_present(const unsigned int type)
{
	struct S_ELEMENT *ELEMENT = DB.ELEMENT;

	while (ELEMENT != NULL) {
		if (type == ELEMENT->type) {
			if (ELEMENT->present)
				return 1;
			else
				return 0;
		}
		ELEMENT = ELEMENT->next;
	}
	printf("Error: Element type not found (present).\n"), exit(1);
}

unsigned int get_Eclass(const unsigned int type)
{
	if (type == POINT || type == LINE || type == QUAD || type == HEX)
		return C_TP;
	else if (type == TRI || type == TET)
		return C_SI;
	else if (type == PYR)
		return C_PYR;
	else if (type == WEDGE)
		return C_WEDGE;

	printf("Error: There is not yet an element class associated with the type provided.\n"), exit(1);
}

struct S_ELEMENT *get_ELEMENT_type(const unsigned int type)
{
	struct S_ELEMENT *ELEMENT = DB.ELEMENT;

	while (ELEMENT != NULL) {
		if (type == ELEMENT->type)
			return ELEMENT;

		ELEMENT = ELEMENT->next;
	}
	printf("Error: Element type not found (type).\n"), exit(1);
}

struct S_ELEMENT *get_ELEMENT_Eclass(const unsigned int type, const unsigned int IndEclass)
{
	struct S_ELEMENT *ELEMENT = DB.ELEMENT;

	if (type == POINT || type == LINE || type == QUAD || type == HEX || (type == WEDGE && IndEclass == 1)) {
		while (ELEMENT != NULL) {
			if (ELEMENT->type == LINE)
				return ELEMENT;

			ELEMENT = ELEMENT->next;
		}
	} else if (type == TRI || type == TET || (type == WEDGE && IndEclass == 0)) {
		while (ELEMENT != NULL) {
			if (ELEMENT->type == TRI)
				return ELEMENT;

			ELEMENT = ELEMENT->next;
		}
	} else if (type == PYR) {
		while (ELEMENT != NULL) {
			if (ELEMENT->type == PYR)
				return ELEMENT;

			ELEMENT = ELEMENT->next;
		}
	}
	printf("Error: Element class not found.\n"), exit(1);
}

unsigned int get_IndFType(const unsigned int Eclass, const unsigned int f)
{
	if (Eclass == C_TP || Eclass == C_SI || (Eclass == C_PYR && f < 4)) {
		return 0;
	} else if (Eclass == C_PYR && f < 5) {
		return 1;
	}
	printf("%d %d\n",Eclass,f);
	printf("Error: Unsupported Eclass/f combination in get_IndFType.\n"), exit(1);
}


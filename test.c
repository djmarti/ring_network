#include "ext_inputs.h"

int main()
{
    /*    char string[100];
       char s[100];
       strcpy(string, "lolailo");
       strcpy(string, "perifollao");
       strcpy(s, "to");
       printf("%s**\n", strcpy(s, string));
       return 0; */

    struct pulse *list;
    list = NULL;
    list = append_pulse(list, new_pulse(0.2, 30, 0.2, 1.2));
    list = append_pulse(list, new_pulse(0.4, 31, 0.1, 1.0));
    print_pulses(list);
    return 0;
}

#ifndef linked_list_h
#define linked_list_h

typedef struct linked_list{
  struct linked_list_member *first, *last;
  int n_members;
} linked_list;

typedef struct linked_list_member {
  struct linked_list_member *prev, *next;
  void * data;
} linked_list_member;

linked_list *
new_linked_list(void);
linked_list *
delete_linked_list(linked_list *head);
linked_list_member *
new_link( linked_list *head );
linked_list_member *
new_link_before( linked_list *head, linked_list_member *old_link );
linked_list_member *
new_last_link( linked_list *head );
void 
delete_link( linked_list_member *this_link, linked_list *head );
void
link_insert_first( linked_list *head, linked_list_member *new_link );

#endif

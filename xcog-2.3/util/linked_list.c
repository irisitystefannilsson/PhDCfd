#include <malloc.h>
#include <stddef.h>
#include <stdlib.h>
#include "linked_list.h"

linked_list *
new_linked_list(void){
  linked_list *head;
  head = (linked_list *) malloc( sizeof(linked_list) );
  head->first = NULL;
  head->last = NULL;
  head->n_members = 0;
  return head;
}

linked_list *
delete_linked_list(linked_list *head){
  linked_list_member *this_link, *next_victim;
  
  if (head == NULL) return NULL;

  this_link = head->first;
  while ( this_link != NULL ){
    next_victim = this_link->next;
/* assume that the structure `data' was pointing to already has been freed */
    free( this_link );
    this_link = next_victim;
  }
  free( head );
  return NULL;
}
  
linked_list_member *
new_link( linked_list *head ){
  linked_list_member *new_link;

/* make a new link */
  new_link = (linked_list_member *) malloc( sizeof(linked_list_member) );
  new_link->prev = NULL;
  new_link->data = NULL;

/* insert the link pointed to by `new_link' at the beginning of the linked list */
  if (head->first != NULL){
    head->first->prev = new_link;
  }
  new_link->next = head->first;
  head->first = new_link;
  if (head->last == NULL)
    head->last = new_link;

/* increase the count of links */
  head->n_members++;

  return new_link;
}

linked_list_member *
new_last_link( linked_list *head ){
  linked_list_member *new_link;

/* make a new link */
  new_link = (linked_list_member *) malloc( sizeof(linked_list_member) );
  new_link->next = NULL;
  new_link->data = NULL;

/* insert the link pointed to by `new_link' at the end of the linked list */
  if (head->last != NULL){
    head->last->next = new_link;
  }
  new_link->prev = head->last;
  head->last = new_link;
  if (head->first == NULL)
    head->first = new_link;

/* increase the count of links */
  head->n_members++;

  return new_link;
}

linked_list_member *
new_link_before( linked_list *head, linked_list_member *old_link ){
  linked_list_member *new_link;

  if (old_link == NULL || head == NULL) return NULL;

/* make a new link */
  new_link = (linked_list_member *) malloc( sizeof(linked_list_member) );
  new_link->prev = NULL;
  new_link->data = NULL;

/* insert the link pointed to by `new_link' before `old_link' */
  new_link->next = old_link;

/* is old_link alone in the existing list? */
  if (head->first == old_link)
    head->first = new_link;
  else{
    old_link->prev->next = new_link;
    new_link->prev = old_link->prev;
  }

  old_link->prev = new_link;

/* increase the count of links */
  head->n_members++;

  return new_link;
}

void 
delete_link( linked_list_member *this_link, linked_list *head ){

/* remove the component grid pointed to by `this_link' from the linked list */
  if (head->first == this_link)
    head->first = this_link->next;
  if (head->last == this_link)
    head->last = this_link->prev;
  if (this_link->prev != NULL)
    this_link->prev->next = this_link->next;
  if (this_link->next != NULL)
    this_link->next->prev = this_link->prev;

/* free the memory occupied by the link */
  free( this_link );

/* increase the count of links */
  head->n_members--;
}

void
link_insert_first( linked_list *head, linked_list_member *new_link ){
  if (head == NULL) return;

/* insert the link pointed to by `new_link' at the beginning of the linked list */
  if (head->first != NULL){
    head->first->prev = new_link;
  }
  new_link->next = head->first;
  head->first = new_link;
  if (head->last == NULL)
    head->last = new_link;

/* increase the count of links */
  head->n_members++;

}


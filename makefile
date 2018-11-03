all:
	gcc wargame_patron.c -o wargame
	gcc alpha-beta.c -o alpha
wargame:
	gcc wargame_patron.c -o wargame
alpha:
	gcc alpha-beta.c -o alpha
clean:
	rm wargame alpha

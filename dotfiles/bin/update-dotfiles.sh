#!/bin/bash
HOME="/home/jy"
DOTFILES="/home/jy/work/src/dotfiles"

mkdir -p "$DOTFILES/conf" "$DOTFILES/misc"

# ~/
cp $HOME/.condarc $DOTFILES/condarc
cp $HOME/.bashrc $DOTFILES/bashrc
cp $HOME/.gitconfig $DOTFILES/gitconfig
cp $HOME/.Rprofile $DOTFILES/Rprofile

# /etc/
cp /etc/pacman.conf $DOTFILES/conf/
cp /etc/fonts/conf.d/64-language-selector-prefer.conf $DOTFILES/conf/

# pacman -Qqen > packages-official.txt 
pacman -Qqen > $DOTFILES/misc/packages-official.txt

# pacman -Qqem > packages-archlinuxcn.txt
pacman -Qqem > $DOTFILES/misc/packages-archlinuxcn.txt
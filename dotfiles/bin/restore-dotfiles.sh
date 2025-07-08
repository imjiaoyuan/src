#!/bin/bash
HOME="/home/jy"
DOTFILES="/home/jy/work/src/dotfiles"

# ~/
cp $DOTFILES/condarc $HOME/.condarc
cp $DOTFILES/bashrc $HOME/.bashrc
cp $DOTFILES/gitconfig $HOME/.gitconfig
cp $DOTFILES/Rprofile $HOME/.Rprofile

# /etc/
cp $DOTFILES/conf/pacman.conf /etc/
cp $DOTFILES/conf/64-language-selector-prefer.conf /etc/fonts/conf.d/

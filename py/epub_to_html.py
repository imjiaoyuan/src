import os
import shutil
import glob
from urllib.parse import urldefrag
import ebooklib
from ebooklib import epub
from bs4 import BeautifulSoup, XMLParsedAsHTMLWarning
import warnings
import re

warnings.filterwarnings("ignore", category=XMLParsedAsHTMLWarning)

INPUT_DIRECTORY = '/home/jy/work/epubs'
MAIN_OUTPUT_DIRECTORY = '/home/jy/work/books'
FAVICON_SOURCE_PATH = '/home/jy/work/books/favicon.ico'
EXCLUDE_FROM_CLEANUP = ['favicon.ico']

def natural_sort_key(filename):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', filename)]

def create_master_index(output_dir, books_list):
    books_html_list = "".join(
        f'<li class="book"><a href="{book["path"]}"><span>{book["title"]}</span></a></li>' for book in books_list
    )
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>My Bookshelf</title>
        <link rel="icon" type="image/x-icon" href="favicon.ico">
        <style>
            body {{
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
                background-color: #ffffff;
                margin: 0;
                padding: 2em;
            }}
            .container {{
                max-width: 1400px;
                margin: 0 auto;
            }}
            h1 {{
                text-align: center;
                font-size: 2.8em;
                margin-bottom: 1.5em;
                color: #2c3e50;
                font-weight: 300;
            }}
            .book-list {{
                display: grid;
                grid-template-columns: repeat(auto-fill, minmax(160px, 1fr));
                gap: 2.5em;
                list-style-type: none;
                padding: 0;
                justify-content: center;
            }}
            .book a {{
                display: flex;
                flex-direction: column;
                justify-content: center;
                align-items: center;
                width: 160px;
                height: 230px;
                margin: 0 auto;
                padding: 1.2em;
                background-color: #ffffff;
                border: 2px solid #e0e0e0;
                border-radius: 5px;
                box-shadow: 0 3px 8px rgba(0, 0, 0, 0.05);
                text-decoration: none;
                color: #333;
                text-align: center;
                font-weight: 600;
                font-size: 1em;
                transition: border-color 0.2s ease, background-color 0.2s ease;
                box-sizing: border-box;
                overflow: hidden;
            }}
            .book a:hover {{
                border-color: #888888;
                background-color: #f9f9f9;
            }}
            .book a span {{
                display: -webkit-box;
                -webkit-box-orient: vertical;
                -webkit-line-clamp: 4;
                line-clamp: 4;
                overflow: hidden;
                word-break: break-word;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>My Bookshelf</h1>
            <ul class="book-list">
                {books_html_list}
            </ul>
        </div>
    </body>
    </html>
    """
    with open(os.path.join(output_dir, 'index.html'), 'w', encoding='utf-8') as f:
        f.write(html_content)

def add_navigation_buttons(soup, prev_chapter, next_chapter, book_title):
    if not soup.body:
        return
    
    nav_style = """
    .chapter-navigation {
        margin-top: 3em;
        padding: 2em 0;
        border-top: 1px solid #ddd;
        text-align: center;
    }
    .chapter-navigation a {
        color: #8e44ad;
        text-decoration: underline;
        margin: 0 1em;
    }
    .chapter-navigation a:hover {
        color: #6b2c91;
    }
    .chapter-navigation .disabled {
        color: #ccc;
        text-decoration: underline;
        margin: 0 1em;
    }
    """
    
    existing_style = soup.find('style')
    if existing_style:
        existing_style.string += nav_style
    else:
        if soup.head:
            style_tag = soup.new_tag('style', string=nav_style)
            soup.head.append(style_tag)
    
    nav_div = soup.new_tag('div', **{'class': 'chapter-navigation'})
    
    if prev_chapter:
        prev_link = soup.new_tag('a', href=prev_chapter)
        prev_link.string = '← 上一章'
        nav_div.append(prev_link)
    else:
        prev_disabled = soup.new_tag('span', **{'class': 'disabled'})
        prev_disabled.string = '← 上一章'
        nav_div.append(prev_disabled)
    
    home_link = soup.new_tag('a', href='../index.html')
    home_link.string = '目录'
    nav_div.append(home_link)
    
    bookshelf_link = soup.new_tag('a', href='../../index.html')
    bookshelf_link.string = '首页'
    nav_div.append(bookshelf_link)
    
    if next_chapter:
        next_link = soup.new_tag('a', href=next_chapter)
        next_link.string = '下一章 →'
        nav_div.append(next_link)
    else:
        next_disabled = soup.new_tag('span', **{'class': 'disabled'})
        next_disabled.string = '下一章 →'
        nav_div.append(next_disabled)
    
    soup.body.append(nav_div)

def convert_ebook_to_html(epub_path, output_dir):
    book = epub.read_epub(epub_path)
    try:
        book_title = book.get_metadata('DC', 'title')[0][0]
    except (IndexError, KeyError):
        book_title = os.path.splitext(os.path.basename(epub_path))[0]

    chapters_dir = os.path.join(output_dir, 'chapters')
    resources_dir = os.path.join(output_dir, 'resources')
    os.makedirs(chapters_dir, exist_ok=True)
    os.makedirs(resources_dir, exist_ok=True)

    chapter_files = []
    for item in book.get_items():
        if item.get_type() == ebooklib.ITEM_DOCUMENT:
            chapter_files.append(os.path.basename(item.get_name()))
    
    chapter_files.sort(key=natural_sort_key)
    chapter_index = {filename: i for i, filename in enumerate(chapter_files)}

    title_map = {}
    def build_title_map(toc_items):
        for item in toc_items:
            link_item, children = (item, None) if isinstance(item, epub.Link) else item
            book_item = book.get_item_with_href(urldefrag(link_item.href).url)
            if book_item:
                title_map[os.path.basename(book_item.get_name())] = link_item.title
            if children: build_title_map(children)
    build_title_map(book.toc)

    for item in book.get_items():
        filename = item.get_name()
        if item.get_type() == ebooklib.ITEM_DOCUMENT:
            content = item.get_content().decode('utf-8', 'ignore')
            soup = BeautifulSoup(content, 'lxml')
            
            page_title = title_map.get(os.path.basename(filename), book_title)
            if soup.title: soup.title.string = page_title
            elif soup.head: soup.head.append(soup.new_tag('title', string=page_title))
                
            if soup.head:
                style_string = """
                body {
                    max-width: 75%; 
                    margin: 2em auto; 
                    line-height: 1.8; 
                    background-color: #e8f5e9; 
                    color: #1b5e20;
                }
                body a {
                    color: #8e44ad;
                }
                """
                style_tag = soup.new_tag('style', string=style_string)
                soup.head.append(style_tag)
                favicon_link = soup.new_tag('link', rel="icon", type="image/x-icon", href="../../favicon.ico")
                soup.head.append(favicon_link)
            
            original_html_dir = os.path.dirname(filename)
            for tag in soup.find_all(['link', 'img']):
                attr = 'href' if tag.has_attr('href') else 'src'
                if tag.has_attr(attr) and not tag[attr].startswith(('http', '#', '/')):
                    original_link = tag[attr]
                    normalized_path = os.path.normpath(os.path.join(original_html_dir, original_link))
                    tag[attr] = f"../resources/{normalized_path.replace(os.sep, '/')}"
            
            current_filename = os.path.basename(filename)
            if current_filename in chapter_index:
                current_index = chapter_index[current_filename]
                prev_chapter = chapter_files[current_index - 1] if current_index > 0 else None
                next_chapter = chapter_files[current_index + 1] if current_index < len(chapter_files) - 1 else None
                
                add_navigation_buttons(soup, prev_chapter, next_chapter, book_title)
            
            with open(os.path.join(chapters_dir, os.path.basename(filename)), 'w', encoding='utf-8') as f:
                f.write(str(soup))
        else:
            res_path = os.path.join(resources_dir, filename)
            os.makedirs(os.path.dirname(res_path), exist_ok=True)
            with open(res_path, 'wb') as f: f.write(item.get_content())
    
    toc_links = []
    def parse_toc(toc_items):
        for item in toc_items:
            link_item, children = (item, None) if isinstance(item, epub.Link) else item
            book_item = book.get_item_with_href(urldefrag(link_item.href).url)
            if book_item and book_item.get_type() == ebooklib.ITEM_DOCUMENT:
                href = f"chapters/{os.path.basename(book_item.get_name())}"
                toc_links.append(f'<li><a href="{href}">{link_item.title}</a></li>')
                if children:
                    toc_links.append('<ul>')
                    parse_toc(children)
                    toc_links.append('</ul>')
                toc_links.append('</li>')
    parse_toc(book.toc)
    
    index_html = f"""
    <!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">
    <title>Table of Contents - {book_title}</title>
    <link rel="icon" type="image/x-icon" href="../favicon.ico">
    <style>
    body {{font-family:sans-serif;max-width:800px;margin:auto;padding:2em;}}
    h1 {{text-align:center;}} ul {{list-style:none;padding-left:1em;}}
    a {{text-decoration:none;color:#0056b3;}} a:hover {{text-decoration:underline;}}
    .toc-header {{display:flex;justify-content:space-between;align-items:center;}}
    </style></head><body><h1>{book_title}</h1><div class="toc-header"><h2>Table of Contents</h2><a href="../index.html">首页</a></div><ul>{"".join(toc_links)}</ul></body></html>
    """
    with open(os.path.join(output_dir, 'index.html'), 'w', encoding='utf-8') as f:
        f.write(index_html)
    return book_title

if __name__ == '__main__':
    if not os.path.isdir(INPUT_DIRECTORY):
        print(f"Error: Input directory does not exist -> {INPUT_DIRECTORY}")
        exit()

    if os.path.isdir(MAIN_OUTPUT_DIRECTORY):
        for item_name in os.listdir(MAIN_OUTPUT_DIRECTORY):
            if item_name not in EXCLUDE_FROM_CLEANUP:
                item_path = os.path.join(MAIN_OUTPUT_DIRECTORY, item_name)
                try:
                    if os.path.isdir(item_path):
                        shutil.rmtree(item_path)
                    else:
                        os.remove(item_path)
                except OSError as e:
                    print(f"Warning: Could not remove {item_path}. Reason: {e}")
    else:
        os.makedirs(MAIN_OUTPUT_DIRECTORY, exist_ok=True)
    
    favicon_dest_path = os.path.join(MAIN_OUTPUT_DIRECTORY, 'favicon.ico')
    if os.path.exists(FAVICON_SOURCE_PATH) and not os.path.exists(favicon_dest_path):
        shutil.copy(FAVICON_SOURCE_PATH, favicon_dest_path)
        
    ebook_files = sorted(glob.glob(os.path.join(INPUT_DIRECTORY, '*.epub')))

    if not ebook_files:
        print(f"No .epub files found in '{INPUT_DIRECTORY}'.")
        exit()

    converted_books = []
    for i, ebook_file_path in enumerate(ebook_files):
        book_name = os.path.basename(ebook_file_path)
        print(f"[{i+1}/{len(ebook_files)}] Converting: {book_name}")
        
        output_book_dir_name = os.path.splitext(book_name)[0]
        final_output_path = os.path.join(MAIN_OUTPUT_DIRECTORY, output_book_dir_name)
        
        try:
            book_title = convert_ebook_to_html(ebook_file_path, final_output_path)
            converted_books.append({'title': book_title, 'path': f"{output_book_dir_name}/index.html"})
            print(f"  -> Success: Output to {final_output_path}")
        except Exception as e:
            print(f"  -> Failed: {e}")
    
    create_master_index(MAIN_OUTPUT_DIRECTORY, converted_books)
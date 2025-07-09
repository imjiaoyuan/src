import os
import shutil
import glob
from urllib.parse import urldefrag
import ebooklib
from ebooklib import epub
from bs4 import BeautifulSoup, XMLParsedAsHTMLWarning
import warnings

warnings.filterwarnings("ignore", category=XMLParsedAsHTMLWarning)

INPUT_DIRECTORY = '/home/jy/work/epubs'
MAIN_OUTPUT_DIRECTORY = '/home/jy/work/books'

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
                style = soup.new_tag('style', string="body {max-width:75%;margin:2em auto;line-height:1.8;}")
                soup.head.append(style)
            
            original_html_dir = os.path.dirname(filename)
            for tag in soup.find_all(['link', 'img']):
                attr = 'href' if tag.has_attr('href') else 'src'
                if tag.has_attr(attr) and not tag[attr].startswith(('http', '#', '/')):
                    original_link = tag[attr]
                    normalized_path = os.path.normpath(os.path.join(original_html_dir, original_link))
                    tag[attr] = f"../resources/{normalized_path}"
            
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
    <title>Table of Contents - {book_title}</title><style>
    body {{font-family:sans-serif;max-width:800px;margin:auto;padding:2em;}}
    h1 {{text-align:center;}} ul {{list-style:none;padding-left:1em;}}
    a {{text-decoration:none;color:#0056b3;}} a:hover {{text-decoration:underline;}}
    </style></head><body><h1>{book_title}</h1><h2>Table of Contents</h2><ul>{"".join(toc_links)}</ul></body></html>
    """
    with open(os.path.join(output_dir, 'index.html'), 'w', encoding='utf-8') as f:
        f.write(index_html)
    return book_title

if __name__ == '__main__':
    if not os.path.isdir(INPUT_DIRECTORY):
        print(f"Error: Input directory does not exist -> {INPUT_DIRECTORY}")
        exit()

    if os.path.exists(MAIN_OUTPUT_DIRECTORY):
        shutil.rmtree(MAIN_OUTPUT_DIRECTORY)
    os.makedirs(MAIN_OUTPUT_DIRECTORY, exist_ok=True)
    
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
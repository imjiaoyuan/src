from markdown import markdown
import jinja2
from pathlib import Path
import re

def read_markdown_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as f:
        return f.read()

def read_template_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as f:
        return f.read()

def convert_markdown_to_html(md_content):
    md_content = re.sub(r'^(\s{2,})(?=\S)', lambda m: '    ' * (len(m.group(1)) // 2), md_content, flags=re.MULTILINE)
    md_content = re.sub(r'<!--.*?-->', '', md_content, flags=re.DOTALL)
    
    extensions = ['extra', 'tables']
    html_content = markdown(md_content, extensions=extensions)
    
    html_content = html_content.replace('<hr>', '<hr style="border: none; border-top: 1px solid #000000; margin: 6px 0;">')
    bold_style = 'style="color: #000000; font-weight: 900; letter-spacing: 0.02em;"'
    html_content = html_content.replace('<strong>', f'<strong {bold_style}>')
    
    html_content = f'<div style="font-size: 0.9em;">{html_content}</div>'
    
    return html_content

def generate_html(md_file_path, output_path, template_path):
    html_content = convert_markdown_to_html(read_markdown_file(md_file_path))
    template = jinja2.Template(read_template_file(template_path))
    html_output = template.render(content=html_content)

    html_path = Path(output_path).with_suffix('.html')
    html_path.parent.mkdir(parents=True, exist_ok=True)
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_output)
    print(f"HTML生成成功: {html_path}")

if __name__ == '__main__':
    cv_dir = '/home/jy/work/src/cv'
    input_file = f'{cv_dir}/cv.md'
    output_file = f'{cv_dir}/cv.html'
    template_file = f'{cv_dir}/template.html'
    
    generate_html(
        md_file_path=input_file,
        output_path=output_file,
        template_path=template_file
    )
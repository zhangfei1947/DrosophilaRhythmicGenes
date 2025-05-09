# app.py
import os, re
import pickle
import json
import datetime
import pytz
from flask import Flask, request, jsonify, render_template, send_from_directory
import lmdb
import csv
import requests

app = Flask(__name__)

# --------------------------
# 加载基因注释
# --------------------------
GENE_ANNOTATION_PATH = "data/gene_annotation.csv"
gene_id_to_name = {}

with open(GENE_ANNOTATION_PATH) as f:
    reader = csv.reader(f)
    next(reader)  # 跳过标题
    for row in reader:
        gene_id_to_name[row[0]] = row[1]

# --------------------------
# 加载4个 LMDB 环境
# --------------------------
SAMPLES = ["T25LD", "T18DD", "T25DD", "T29DD"]
lmdb_envs = {
    sample: lmdb.open(os.path.join("lmdbs", f"{sample}.lmdb"), readonly=True)
    for sample in SAMPLES
}


# --------------------------
# 访客记录 (SQLite)
# --------------------------
def get_central_time():
    central = pytz.timezone('America/Chicago')  # 中部时间，自动处理夏令时
    return datetime.datetime.now(central)

from flask_sqlalchemy import SQLAlchemy
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///visitors.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

class Visitor(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    ip = db.Column(db.String(15), nullable=False)
    country = db.Column(db.String(100))
    region = db.Column(db.String(100))
    visit_time = db.Column(db.DateTime, default=get_central_time)
    query = db.Column(db.String(500))

@app.before_request
def track_visitor():
    # 只记录 /gene 路由的请求
    if request.endpoint == 'get_gene_expr':
        # 获取真实 IP：优先使用 X-Forwarded-For
        x_forwarded_for = request.headers.get('X-Forwarded-For')
        if x_forwarded_for:
            ip = x_forwarded_for.split(',')[0].strip()  # 第一个 IP 是客户端的真实 IP
        else:
            ip = request.remote_addr  

        query = request.args.get('gene_ids', '')
        
        # 获取国家（带异常处理）
        country = "Unknown"
        try:
            res = requests.get(f'https://ipapi.co/{ip}/json/', timeout=3).json()
            country = res.get('country_name', 'Unknown')
            region = res.get('region', 'Unknown')
        except Exception as e:
            print("Error fetching country:", e)
            country = 'Unknown'
            region = 'Unknown'

        # 存入数据库
        try:
            visitor = Visitor(ip=ip, country=country, region=region, query=query)
            db.session.add(visitor)
            db.session.commit()
        except Exception as e:
            db.session.rollback()
            print("Database error:", e)

# --------------------------
# 请求频率限制
# --------------------------
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

limiter = Limiter(
    app=app,
    key_func=get_remote_address,
    default_limits=["200 per day", "10 per minute"]
)

# --------------------------
# JSON API 接口
# --------------------------
@app.route('/gene')
@limiter.limit("10 per minute")
def get_gene_expr():
    gene_ids = re.sub(r'\s+', '', request.args.get('gene_ids', '')).split(',')

    gene_ids = [g.strip() for g in gene_ids if g.strip()]

    if not gene_ids:
        return jsonify({"error": "No gene IDs provided"}), 400
    
    if len(gene_ids) > 100:
        return jsonify({"error": "Maximum 100 gene IDs allowed"}), 400

    results = []

    for gene_id in gene_ids:
        if gene_id not in gene_id_to_name:
            continue

        gene_data = {"gene_id": gene_id, "gene_name": gene_id_to_name[gene_id], "expression": {}}

        for sample in SAMPLES:
            with lmdb_envs[sample].begin() as txn:
                data = txn.get(gene_id.encode())
                if data:
                    fpkms = pickle.loads(data)
                    times = sorted(fpkms.keys())
                    gene_data["expression"][sample] = {
                        "time_points": times,
                        "rep1": [fpkms[t]["r1"] for t in times],
                        "rep2": [fpkms[t]["r2"] for t in times],
                        "mean": [fpkms[t]["mean"] for t in times]
                    }

        if gene_data["expression"]:
            results.append(gene_data)

    return jsonify(results)

# --------------------------
# 页面路由
# --------------------------
@app.route('/')
def index():
    return render_template("index.html")

# --------------------------
# 文件下载
# --------------------------
@app.route('/download')
def download():
    return send_from_directory("data", "GeneExpressionData.zip", as_attachment=True)

# --------------------------
# 访客记录页面
# --------------------------
@app.route('/visitors')
def visitors():
    page = request.args.get('page', 1, type=int)
    per_page = 20  # 每页显示20条记录
    
    # 使用新的查询语法
    stmt = db.select(Visitor).order_by(Visitor.visit_time.desc())
    pagination = db.paginate(stmt, page=page, per_page=per_page, error_out=False)
    visitors_list = pagination.items
    
    # 获取统计数据
    total_visits = db.session.scalar(db.select(db.func.count()).select_from(Visitor))
    unique_ips = db.session.scalar(db.select(db.func.count(db.distinct(Visitor.ip))))
    
    # 获取国家统计
    stmt = db.select(Visitor.country, db.func.count(Visitor.id)).\
        group_by(Visitor.country).\
        order_by(db.func.count(Visitor.id).desc())
    countries = db.session.execute(stmt).all()
    
    # 获取最近24小时的访问量
    yesterday = datetime.datetime.utcnow() - datetime.timedelta(days=1)
    stmt = db.select(db.func.count()).select_from(Visitor).where(Visitor.visit_time >= yesterday)
    recent_visits = db.session.scalar(stmt)
    
    # 获取最常查询的基因
    top_genes = []
    if total_visits > 0:
        gene_counts = {}
        for visitor in db.session.execute(db.select(Visitor)).scalars():
            if visitor.query:
                genes = visitor.query.split(',')
                for gene in genes:
                    gene = gene.strip()
                    if gene:
                        gene_counts[gene] = gene_counts.get(gene, 0) + 1
        
        top_genes = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)[:10]
    
    return render_template('visitors.html', 
                          visitors=visitors_list,
                          pagination=pagination,
                          total_visits=total_visits,
                          unique_ips=unique_ips,
                          countries=countries,
                          recent_visits=recent_visits,
                          top_genes=top_genes)


# --------------------------
# 初始化数据库
# --------------------------
with app.app_context():
    db.create_all()

if __name__ == '__main__':
    app.run(debug=True)

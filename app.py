# app.py
import os
import pickle
import json
import datetime
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
from flask_sqlalchemy import SQLAlchemy
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///visitors.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

class Visitor(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    ip = db.Column(db.String(15), nullable=False)
    country = db.Column(db.String(100))
    visit_time = db.Column(db.DateTime, default=datetime.datetime.utcnow)
    query = db.Column(db.String(500))

@app.before_request
def track_visitor():
    if request.endpoint == 'get_gene_expr':
        ip = request.remote_addr
        query = request.args.get('gene_ids', '')
        
        try:
            res = requests.get(f'https://ipapi.co/{ip}/json/', timeout=2).json()
            country = res.get('country_name', 'Unknown')
        except:
            country = 'Unknown'

        visitor = Visitor(ip=ip, country=country, query=query)
        db.session.add(visitor)
        db.session.commit()

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
    gene_ids = request.args.get('gene_ids', '').split(',')
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
    
    # 获取分页数据
    pagination = Visitor.query.order_by(Visitor.visit_time.desc()).paginate(
        page=page, per_page=per_page, error_out=False)
    visitors_list = pagination.items
    
    # 获取统计数据
    total_visits = Visitor.query.count()
    unique_ips = db.session.query(Visitor.ip).distinct().count()
    countries = db.session.query(Visitor.country, db.func.count(Visitor.id)).\
        group_by(Visitor.country).order_by(db.func.count(Visitor.id).desc()).all()
    
    # 获取最近24小时的访问量
    yesterday = datetime.datetime.utcnow() - datetime.timedelta(days=1)
    recent_visits = Visitor.query.filter(Visitor.visit_time >= yesterday).count()
    
    # 获取最常查询的基因
    top_genes = []
    if total_visits > 0:
        gene_counts = {}
        for visitor in Visitor.query.all():
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
